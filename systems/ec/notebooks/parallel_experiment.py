from mcm.tcg import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool
import multiprocessing as mp
import importlib
from from_perplex import *

yr = 3.154e7
kyr = 1e3*yr
Myr = 1e6*yr
s = 1
mm = 1e-3
km = 1e3
g = 1e-3
cm = 1e-2

### ------------ INPUTS -------------------
reference= 'parallel_experiment'
composition = 'hacker_2015_md_xenolith'
rxn_name = 'eclogitization_agu17_stx21_rx'

# only phases greater than this fraction will be plotted
phasetol = 1.e-5 # 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9
max_steps = 3e4

# Set initial and final crustal thickness (m)
moho_i = 30.*km # typical crust
moho_f = 70.*km # Altiplano, Tibet

# Calc descent rate of Moho
shortening_rate = 3.0 *mm/yr # m/s

L0 = 100.*km # lithospheric thickness - m
z0 = moho_i # 
descent_rate = z0/L0 * shortening_rate # m/s

# choose a time step
dt = 100.*kyr # sec
max_t = 100.*Myr
times = np.arange(0,max_t+dt,dt)

depth_m = z0 + times*descent_rate
depth_m[depth_m > 70.*km] = 70.* km

# Calc pressure, temperature from depths
# TODO: use a steady-state geotherm
P_range = 0.027 * depth_m/1e3 # assumes 0.027 Gpa/km
T_moho_i = 800. + 273.15
T_moho_f = 850. + 273.15
T_range = T_moho_i + (depth_m - depth_m[0])/(depth_m[-1] - depth_m[0])*(T_moho_f - T_moho_i)
print(P_range)
show_equilibrium = False
processes =  mp.cpu_count()-1

# ------------------------------------------

# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-r", "--rxn_name")

    args = parser.parse_args()

    if args.composition is not None:
        composition = args.composition
    if args.rxn_name is not None:
        rxn_name = args.rxn_name

#====================================================

rxn = get_reaction(rxn_name)
phase_names, endmember_names = get_names(rxn)

rxn.set_parameter('T0', T_range[0]) # K; default 2000 K
mi0, Xik0, phii0, Cik0 = get_point_composition(rxn, composition)

ipyrolite = get_rho_interpolator('xu_2008_pyrolite')
iharzburgite = get_rho_interpolator('xu_2008_harzburgite')
rho_pyrolite=ipyrolite((T_range, P_range))
rho_harzburgite=iharzburgite((T_range,P_range))

_Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
_mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0

I = len(rxn.phases())
Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(Kis)

# Equilibrate model at initial (T0, P0)
ode = ScipyPDReactiveODE(rxn)
ode.solve(T_range[0],GPa2Bar(P_range[0]),_mi0,_Cik0,1,Da=1e6,eps=eps)
rho0 = ode.final_rho()*100 # kg/m3

# Damkoehler number
# Da = t0*Gamma0/rho0
# Gamma0 = Da*rho0/t0
# Gamma0 = r0/dg <- grain size
# t0 = dt
Da = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]

# NOTE
# reaction rates are multiplied by Damkholer number
# this implies that they are scaled so that mass flux rate = 1
# this defines a time scale?

# Reaction rates based on Hetenyi et al. 2017
# these are based on available surface area of 3/r where r is grain radius
#reaction_rate_gcm =  [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3] # g/cm2/yr

reaction_rate = [da*rho0/dt for da in Da] # Gamma0 (kg/m3/s)
dg = 1.*mm # grain size
reaction_rate_per_surface = [r*dg for r in reaction_rate] # r0 (kg/m2/s), assumes 1mm grains
reaction_rate_per_surface_gcm = [r/10*yr for r in reaction_rate_per_surface] # g/cm2/yr

print("rho0 = {}".format(rho0))
print("Da = {}".format(Da))
print("rxn rate = {} g/cm2/yr".format(reaction_rate_per_surface_gcm))
print("dt = {} kyr".format(dt/kyr))
print("total t = {} Myr".format(max_t/Myr))

# function to run in parallel
def run_experiment(Da):
    rho_final = np.zeros(T_range.shape)
    phases_final = ['' for i,_ in enumerate(depth_m)]
    phi_final = np.empty(T_range.shape+(I,))
    mi_final = np.zeros(T_range.shape+(I,))
    Cik_final = np.empty(T_range.shape+(K,))
    Xik_final = np.empty(T_range.shape+(K,))

    # Equilibrate model at initial (T0, P0)
    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T_range[0],GPa2Bar(P_range[0]),_mi0,_Cik0,1,Da=1e5,eps=eps)

    # get "equilibrated" composition
    mi1 = ode.sol.y[:ode.I, -1]
    Cik1 = ode.sol.y[ode.I:ode.I+ode.K, -1]

    # Thicken the crust
    for i,t in enumerate(times):
        P = P_range[i]
        T = T_range[i]
        ode = ScipyPDReactiveODE(rxn)

        ode.solve(T, GPa2Bar(P), mi1, Cik1, dt/dt, Da=Da, eps=eps, method="BDF_mcm", max_steps=max_steps)
        odephasenames, phaseabbrev = ode.final_phases(phasetol)
        phases = '+'.join(phaseabbrev)
        rho = ode.final_rho()

        Cik = ode.sol.y[ode.I:ode.I+ode.K,-1]
        mi = ode.sol.y[:ode.I,-1] # -1 = final time step
        mi_reg = ode.regularizem(mi)
        Cik1 = Cik
        mi1 = mi
        C = ode.reshapeC(Cik)
        Xik = rxn.C_to_X(C)
        Cs = ode.regularizeC(C)
        Cik_reg = [c for arr in Cs for c in arr]
        rhoi = rxn.rho(T, GPa2Bar(P), Cs)
        v = mi/rhoi
        vi  = 1./v.sum()
        phi = vi*mi/rhoi
        rho_final[i]=rho
        phases_final[i]=phases 
        phi_final[i]=phi 
        mi_final[i]=mi_reg
        Cik_final[i]=Cik_reg
        Xik_final[i]= [x for arr in Xik for x in arr]

    return rho_final, phases_final, phi_final, mi_final, Cik_final, Xik_final, Da

# Fn to run an "equilibrium" profile
_Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
_mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0
def equilibrate(i):
    P = P_range[i]
    T = T_range[i]
    ode = ScipyPDReactiveODE(rxn)

    ode.solve(T,GPa2Bar(P),_mi0,_Cik0,1e5,Da=1,eps=eps,method="BDF_mcm",max_steps=max_steps)
    odephasenames, phaseabbrev = ode.final_phases(phasetol)
    phases = '+'.join(phaseabbrev)
    rho = ode.final_rho()

    Cik = ode.sol.y[ode.I:ode.I+ode.K,-1]
    mi = ode.sol.y[:ode.I,-1] # -1 = final time step

    C = ode.reshapeC(Cik)
    Cs = ode.regularizeC(C)
    Cik_reg = [c for arr in Cs for c in arr]
    rhoi = rxn.rho(T, GPa2Bar(P), Cs)
    v = mi/rhoi
    vi  = 1./v.sum()
    phi = vi*mi/rhoi

    return rho, phases, phi, mi, Cik_reg, i

# run for varying damkhoeler numbers
with Pool(processes) as pool:
    sols = pool.map(run_experiment, Da)

outputPath = Path("figs",reference,composition,rxn_name)
outputPath.mkdir(parents=True, exist_ok=True)

depth_km = depth_m/1e3

# plot 
num_subplots = 2 + len(phase_names) + 2
fig = plt.figure(figsize=(3*num_subplots,12))
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
ax = fig.add_subplot(1,num_subplots,(1,2))
ax.invert_yaxis()

n = len(Da)
colormap = plt.cm.get_cmap('Greys')
ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))
for rho, *rest in sols:
    plt.plot(rho/10, depth_km)

# Optionally run for equilibrium
if show_equilibrium:
    with Pool(processes,maxtasksperchild=12) as pool:
        sols_eq = pool.map(equilibrate, range(len(P_range)))

    # plot
    rho_equil = np.zeros(T_range.shape)  
    for rho, phases, phi, mi, Cik, Xik, i in sols_eq:
        rho_equil[i]=rho

    plt.plot(rho_equil/10,depth_km,'k--')


plt.plot(rho_pyrolite/10,depth_km, 'r:')
plt.plot(rho_harzburgite/10,depth_km, 'm:')

plt.legend(['$r_0 = ${:.1e} g/cm$^2$/yr'.format(r) for r in reaction_rate_per_surface_gcm] + (['equil.', 'pyrolite', 'harzburgite'] if show_equilibrium else  ['pyrolite', 'harzburgite']), loc="upper right")

plt.ylabel('Depth (km)')
plt.xlabel('Density')
plt.xlim([2.8, 3.8])

# plot all phases

cmaps = [plt.cm.get_cmap(name) for name in ['Blues', 'YlOrBr', 'Greens','Reds','Purples','copper_r']]
for i,phase in enumerate(phase_names):
    ax = fig.add_subplot(1,num_subplots,3+i)
    ax.invert_yaxis()
    cmap = cmaps[i]
    ax.set_prop_cycle(plt.cycler('color', cmap(np.linspace(0.2, 1, n))))
    ax.set_xlim([0., 80.])
    ax.set_ylabel(None)
    ax.set_xlabel("{} (wt%)".format(phase))
    ax.set_xticks(np.arange(0,110,10))
    ax.set_xticklabels([None, None, 20, None, 40, None, 60, None, 80, None, None])
    ax.set_yticklabels([])
    for rho, phases, phi, mi, Cik, Xik, *rest in sols:
        # plot garnet mass fraction
        plt.plot(mi[:,i]*100, depth_km)
    
# Plot plagioclase composition
ax = fig.add_subplot(1,num_subplots,num_subplots-1)
ax.invert_yaxis()
ax.set_xlim([0., 1.0])
ax.set_ylabel(None)
ax.set_yticklabels([])
ax.set_xlabel("$X_{\mathrm{An}}$")
colormap = plt.cm.get_cmap('Greys')
ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))

for rho, phases, phi, mi, Cik, Xik, *rest in sols:
    Xan = np.asarray(Xik[:,10])
    plt.plot(Xan, depth_km)

# Plot Cpx composition
ax = fig.add_subplot(1,num_subplots,num_subplots)
ax.invert_yaxis()
ax.set_xlim([0., 1.0])
ax.set_ylabel(None)
ax.set_yticklabels([])
ax.set_xlabel("$X_{\mathrm{Jd}}$")
colormap = plt.cm.get_cmap('Greys')
ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))

for rho, phases, phi, mi, Cik, Xik, *rest in sols:
    Xjd = np.asarray(Xik[:,4])
    plt.plot(Xjd, depth_km)

fig.suptitle('$R_m=${:.1f} km/Myr, $T_m=${:.0f}-{:.0f} °C, $d_m=${:.0f}-{:.0f} km, {}'.format(descent_rate*1e3,T_moho_i-273.15,T_moho_f-273.15,moho_i/1e3,moho_f/1e3,composition),y=0.9)
plt.savefig(Path(outputPath,"{}.{}".format("results", "pdf")))
plt.savefig(Path(outputPath,"{}.{}".format("results", "png")))
