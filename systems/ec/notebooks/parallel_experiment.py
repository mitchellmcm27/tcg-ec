from mcm import EcModel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool
import multiprocessing as mp
import importlib
import from_perplex as pp

yr = 1 # 3.154e7
mm = 1e-3
g = 1e-3
cm = 1e-2

### ------------ INPUTS -------------------
reference= 'parallel_experiment'
composition = 'hacker_2015_md_xenolith'
rxn_name = 'eclogitization_agu17_stx21_rx'

# only phases greater than this fraction will be plotted
phasetol = 1.e-2 # 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9
max_steps = 1e4

# Set initial and final crustal thickness (km)
moho_i = 30. # typical crust
moho_f = 70. # Altiplano, Tibet

# Calc descent rate of Moho
shortening_rate_kmMyr = 3.0 # km/Myr = mm/yr
shortening_rate_myr = shortening_rate_kmMyr/1e3 # m/yr

L0 = 100. # lithospheric thickness - km

z0 = moho_i # km
descent_rate_kmMyr = z0/L0 * shortening_rate_kmMyr # mm/yr
descent_rate = descent_rate_kmMyr*mm/yr  # = m/yr
avg_density = 2800 # kg/m3, assumes a linear increase in density from 2600 at the surface to 3000 at the moho

# choose a time step
dt = 100.e3 # yrs
depth_step_size = dt*descent_rate # meters

depth_m = np.arange(moho_i*1e3, moho_f*1e3 + depth_step_size, depth_step_size) # meters
depth_km = depth_m/1e3

# Calc pressure, temperature from depths
# TODO: use a steady-state geotherm
P_range = 0.027 * depth_km # GPa
T_moho_i = 800.
T_moho_f = 800.
T_range = T_moho_i + (depth_km - depth_km[0])/(depth_km[-1] - depth_km[0])*(T_moho_f - T_moho_i) + 273. 

# Damkoehler number
Da = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]

# Reaction rates based on Hetenyi et al. 2017
# these are based on available surface area of 3/r where r is grain radius
#reaction_rate_gcm =  [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3] # g/cm2/yr

S0 = 6000 # m2/m3; 1 mm cubes 
h0 = 1e-3 # 1mm
                                #kg/m3      m/yr       
reaction_rate_per_surface = [da*avg_density*shortening_rate_myr for da in Da]# kg/m2/yr
reaction_rate_per_surface_gcm = [r/10 for r in reaction_rate_per_surface] # g/cm2/yr

print("Da = {}".format(Da))
print("rxn rate = {} g/cm2/yr".format(reaction_rate_per_surface_gcm))
print("dt = {} kyr".format(dt))
print("total t = {} Myr".format(dt * len(depth_km)/1e6))

show_equilibrium = False
processes =  mp.cpu_count()

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


rxn = EcModel.get_reaction(rxn_name)

mi0, Xik0, phii0, Cik0 = pp.get_point_composition(composition)

def x2c(rxn, Xik0):
    return np.asarray([c for (i, ph) in enumerate(rxn.phases()) for c in ph.x_to_c(Xik0[i])])

def phi2m(rxn, phii0, Cik0, T=900.,p=10000.):
    '''Converts phase modes in volume fraction to mass fraction given an intial EM composition in mass fractions.'''    
    densities = []
    C = rxn.zero_C()
    Ki = 0
    for i,ph in enumerate(rxn.phases()):
        n = len(ph.endmembers())
        C[i] = Cik0[Ki:Ki+n]
        Ki = Ki+n

    C = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    C = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    densities = [ph.rho(T, p, C[i]) for i,ph in enumerate(rxn.phases())]
    mass = np.sum(np.asarray(densities) * np.asarray(phii0))
    mi0 = np.asarray([v*densities[i]/mass for (i, v) in enumerate(phii0)])

    return mi0

# function to run in parallel
def run_experiment(Da):

    _Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
    _mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0

    I = len(rxn.phases())
    Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
    K = sum(Kis)

    rho_final = np.zeros(T_range.shape)
    s_final = np.zeros(T_range.shape)
    phases_final = ['' for i,_ in enumerate(depth_km)]
    phi_final = np.empty(T_range.shape+(I,))
    mi_final = np.zeros(T_range.shape+(I,))
    Cik_final = np.empty(T_range.shape+(K,))
    Xik_final = np.empty(T_range.shape+(K,))

    # Equilibrate model at initial (T0, P0)
    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T_range[0],GPa2Bar(P_range[0]),_mi0,_Cik0,1e5,Da=1e4,eps=eps)

    # get "equilibrated" composition
    mi1 = ode.sol.y[:ode.I, -1]
    Cik1 = ode.sol.y[ode.I:ode.I+ode.K, -1]

    # Thicken the crust
    for i,P in enumerate(P_range):
        T = T_range[i]
        ode = ScipyPDReactiveODE(rxn)

        ode.solve(T,GPa2Bar(P),mi1,Cik1,dt,Da=Da,eps=eps,method="BDF_mcm",max_steps=max_steps)
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
        s_final[i] = np.sum(rxn.s(T,GPa2Bar(P),C))

    return rho_final, phases_final, phi_final, mi_final, Cik_final, Xik_final, s_final, Da

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

# plot 
num_subplots = 2 + len(phase_names) + 3
fig = plt.figure(figsize=(3*num_subplots,12))
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
ax = fig.add_subplot(1,num_subplots,(1,2))
ax.invert_yaxis()

n = len(Da)
colormap = plt.cm.get_cmap('Greys')
ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))
for rho, phases, phi, mi, Cik, Xik, s, _Da in sols:
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

ipyrolite = pp.get_rho_interpolator('xu_2008_pyrolite')
rho_pyrolite=ipyrolite((T_range, GPa2Bar(P_range)))
plt.plot(rho_pyrolite/10,depth_km, 'r:')

iharzburgite = pp.get_rho_interpolator('xu_2008_harzburgite')
rho_harzburgite=iharzburgite((T_range,GPa2Bar(P_range)))
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

# Plot Specific entropy
#ax = fig.add_subplot(1,num_subplots,num_subplots)
#ax.invert_yaxis()
#ax.set_ylabel(None)
#ax.set_xlabel("specific $S$")
#
#colormap = plt.cm.get_cmap('Greens')
#ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))
#
#for rho, phases, phi, mi, Cik, Xik, s, _Da in sols:
#    # s is J/K/kg
#    plt.plot(s, depth_km)

fig.suptitle('$R_m=${:.1f} km/Myr, $T_m=${:.0f}-{:.0f} °C, $d_m=${:.0f}-{:.0f} km, {}'.format(descent_rate_kmMyr,T_moho_i, T_moho_f, moho_i,moho_f,composition),y=0.9)
plt.savefig(Path(outputPath,"{}.{}".format("results", "pdf")))
