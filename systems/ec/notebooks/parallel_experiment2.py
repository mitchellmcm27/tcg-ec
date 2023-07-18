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
from scipy.integrate import solve_ivp
from get_geotherm import *
from geotherm_steady import * 

yr = 3.154e7
kyr = 1e3*yr
Myr = 1e6*yr
s = 1
mm = 1e-3
km = 1e3
g = 1e-3
cm = 1e-2

### ------------ INPUTS -------------------
reference= 'parallel_experiment2'
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

L0 = 60.*km # lithospheric thickness - m
z0 = moho_i # 
descent_rate = z0/L0 * shortening_rate # m/s

# choose a time step
dt = 100.*kyr # sec
max_t = (moho_f - moho_i) / descent_rate # seconds
# Calc pressure, temperature from depths
# TODO: use a steady-state geotherm
T_moho_i = get_geotherm(L0,z0,1,3e7*1e6)
T_moho_f = get_geotherm(L0,z0,moho_f/moho_i,3e7*1e6)
print(T_moho_f)
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

# initial temperature and pressure in K and bars
T0 = T_moho_i
P0 = 2780 * 9.81 * 30.e3/1e5 # Gpa, 0.027 Gpa/km
print(P0)
rxn.set_parameter('T0', T0) # K; default 2000 K
mi0, Xik0, phii0, Cik0 = get_point_composition(rxn, composition)

ipyrolite = get_rho_interpolator('xu_2008_pyrolite')
iharzburgite = get_rho_interpolator('xu_2008_harzburgite')

_Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
_mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0

I = len(rxn.phases())
Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(Kis)

# Equilibrate model at initial (T0, P0)
ode = ScipyPDReactiveODE(rxn)
ode.solve(T0, P0,_mi0,_Cik0,1,Da=1e6,eps=eps)
rho0 = ode.final_rho()*100 # kg/m3
Cik0 = ode.sol.y[ode.I:ode.I+ode.K,-1]
mi0 = ode.sol.y[:ode.I,-1] # -1 = final time step

# Damkoehler number
Das = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5]

r0 = [da*rho0/max_t for da in Das] # Gamma0 (kg/m3/s)
dg = 3.*mm # grain size
reaction_rate_per_surface = [r*dg for r in r0] # r0 (kg/m2/s), assumes 1mm grains
reaction_rate_per_surface_gcm = [r/10*yr for r in reaction_rate_per_surface] # g/cm2/yr
print(reaction_rate_per_surface_gcm)
_Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase


scale= {'T':T0, 'P':P0, 'rho':rho0, 'h':(moho_f-moho_i), 't':max_t}
print(scale)

def rhs(t,u,rxn,scale,Da,verbose=False):

    # Extract variables
    mi = u[:I]
    Cik = u[I:I+K]
    T = u[I+K] # 1 to 1.something
    P = u[I+K+1] # 1 to 2.something
 
    # scale Temperature and pressure
    Ts = scale['T']*T
    Ps = scale['P']*P

    C = rxn.zero_C()
    # reshape C
    Kis = 0
    for i,Ki in enumerate(_Kis):
        C[i] = Cik[Kis:Kis+Ki]
        Kis = Kis+Ki
    # regularize C
    Cs = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    Cs = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    # phase properties
    rhoi = np.array(rxn.rho(Ts, Ps, Cs))

    #alpha = np.array(rxn.alpha(Ts, Ps, C))
    #cp = np.array(rxn.Cp(Ts, Ps, C))
    #s = np.array(rxn.s(Ts, Ps, C))
    v = np.sum(mi/rhoi)

    # mean properties
    #rho_bar = 1./v
    #alpha_bar = v.dot(alpha)*rho_bar
    #cp_bar = mi.dot(cp)
    #s_bar = mi.dot(s)

    #A = scale['P']/scale['rho']
    
    # regularize m
    mis = np.asarray(mi)
    mis = mi + eps

    Gammai = np.asarray(rxn.Gamma_i(Ts,Ps,Cs,mi))
    gamma_ik = rxn.Gamma_ik(Ts,Ps,Cs,mi)
    Gammaik = np.zeros(K)
    sKi = 0
    for i in range(I):
        for k in range(_Kis[i]):
            Gammaik[sKi+k] = gamma_ik[i][k]
        sKi += _Kis[i]
    
    dudt = np.zeros(u.shape)
    sKi = 0
    for i in range(I):
        dudt[i] = Da*rho0*Gammai[i]*v
        for k in range(_Kis[i]):
            GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
            dudt[I+sKi+k] = Da*rho0*GikcGi*v/mis[i]
        sKi += _Kis[i]
    
  
    # linear temperature
    dTdt_linear = (T_moho_f - T_moho_i)/scale['T']

    dzdt = scale['h'] # m per unit time
    dPdt = 2780*9.81/1e5 * dzdt # bar per unit time
    dPdt = dPdt/scale['P']

    shortening = 1 + 1.3*t
    T_steady, dTdz = geotherm_steady(z0/L0, L0*shortening, shortening)

    dTdt = (T_steady-Ts)*dzdt/scale['T']
    dudt[I+K:] = np.array([dTdt, dPdt])
    return dudt

# function to run in parallel
def run_experiment(Da):
    rxn = get_reaction(rxn_name)
    rxn.set_parameter('T0',T0)

    u0=np.empty(I+K+2)
    u0[:I] = mi0
    u0[I:I+K] = Cik0

    u0[I+K:] = np.array([1., 1.]) # T/T0, P/P0
    
    #Pmax = 2780 * 9.81 * moho_f/1e5 # bar
    #print(Pmax)
    #event = lambda t,y,rxn,scale,Da: Pmax - y[-1]*scale['P'] # when pressure = Pmax
    #event.terminal = True
    sol = solve_ivp(rhs,[0,1],u0,args=(rxn,scale,Da),dense_output=True,method='BDF',rtol=rtol,atol=atol,events=None)
    
    T = sol.y[-2]*scale['T']
    P = sol.y[-1]*scale['P']
    print('{} P_end = {:.2f} Gpa. T_end = {:.2f} K. Used {:n} steps: '.format(sol.message,P[-1]/1e4,T[-1],len(sol.t)))

    return sol, Da

# run for varying damkhoeler numbers
with Pool(processes) as pool:
    # blocks until all finished
    sols = pool.map(run_experiment, Das)

outputPath = Path("figs",reference,composition,rxn_name)
outputPath.mkdir(parents=True, exist_ok=True)


# plot 
num_subplots = 3 + len(phase_names) + 2
subplot_mosaic = [['rho','rho'] + phase_names + ['An', 'Jd', 'T']]

fig = plt.figure(figsize=(3*num_subplots,12))
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axes = fig.subplot_mosaic(subplot_mosaic)
[ax.invert_yaxis() for label,ax in axes.items()]

n = len(Das)
z = np.linspace(0,1,100)
depth_m = moho_i + z*(moho_f-moho_i)
depth_km = depth_m / 1000.
greys = plt.cm.get_cmap('Greys')
axes['rho'].set_prop_cycle(plt.cycler('color', greys(np.linspace(0.2, 1, n))))
axes['An'].set_prop_cycle(plt.cycler('color', greys(np.linspace(0.2, 1, n))))
axes['Jd'].set_prop_cycle(plt.cycler('color', greys(np.linspace(0.2, 1, n))))
for i,obj in enumerate(sols):
    ax = axes['rho']

    sol, Da = obj
    y = sol.sol(z)
    mi_times  = y[:I].T
    Cik_times = y[I:I+K].T
    T = y[I+K].T*scale['T']
    P = y[I+K+1].T*scale['P']

    C_times = [ode.reshapeC(Cik) for Cik in Cik_times]
    Cs_times = [ode.regularizeC(C) for C in C_times]
    rho = [1/sum(mi_times[t]/rxn.rho(T[t], P[t], Cs))/10 for t,Cs in enumerate(Cs_times)]
    ax.plot(rho, depth_km)

rho_pyrolite=ipyrolite((T, P/1e4))
rho_harzburgite=iharzburgite((T,P/1e4))

ax.plot(rho_pyrolite/10,depth_km, 'r:')
ax.plot(rho_harzburgite/10,depth_km, 'm:')
ax.legend(['$r_0 = ${:.1e} g/cm$^2$/yr'.format(r) for r in reaction_rate_per_surface_gcm] +  (['pyrolite', 'harzburgite']), loc="upper right")

ax.set_ylabel('Depth (km)')
ax.set_xlabel('Density')
ax.set_xlim([2.8, 3.8])

# Plot temperature profile
ax = axes['T']
ax.plot(T-273.15, depth_km, linewidth=2)
ax.set_xlabel("T (°C)")
ax2 = ax.twinx()
ax2.plot(T-273.15, P/1e4, alpha=0)
ax2.invert_yaxis()
ax2.set_ylabel("P (GPa)")

# plot all phases
cmaps = [plt.cm.get_cmap(name) for name in ['Blues', 'YlOrBr', 'Greens','Reds','Purples','copper_r']]
for i,phase in enumerate(phase_names):
    ax = axes[phase]
    cmap = cmaps[i]
    ax.set_prop_cycle(plt.cycler('color', cmap(np.linspace(0.2, 1, n))))
    ax.set_xlim([0., 80.])
    ax.set_ylabel(None)
    ax.set_xlabel("{} (wt%)".format(phase))
    ax.set_xticks(np.arange(0,110,10))
    ax.set_xticklabels([None, None, 20, None, 40, None, 60, None, 80, None, None])
    ax.set_yticklabels([])
    for j, obj in enumerate(sols):
        sol, Da = obj
        y = sol.sol(z)
        mi_times  = y[i].T
        Cik_times = y[I:I+K].T
        T = y[I+K].T*scale['T']
        P = y[I+K+1].T*scale['P']
        # plot garnet mass fraction
        ax.plot(mi_times*100, depth_km)

# Plot plagioclase composition
ax = axes['An']
ax.set_xlim([0., 1.0])
ax.set_ylabel(None)
ax.set_yticklabels([])
ax.set_xlabel("$X_{\mathrm{An}}$")

ax = axes['Jd']
ax.set_xlim([0., 1.0])
ax.set_ylabel(None)
ax.set_yticklabels([])
ax.set_xlabel("$X_{\mathrm{Jd}}$")


for j, obj in enumerate(sols):
    sol, Da = obj
    y = sol.sol(z)
    mi_times  = y[i].T
    Cik_times = y[I:I+K].T
    C_times = [ode.reshapeC(Cik) for Cik in Cik_times]
    Cs_times = [ode.regularizeC(C) for C in C_times]
    T = y[I+K].T*scale['T']
    P = y[I+K+1].T*scale['P']
    Xik_nested_times = np.asarray([rxn.C_to_X(c) for c in Cs_times])
    XAn = [x[3][0] for x in Xik_nested_times]
    XJd = [x[0][4] for x in Xik_nested_times]

    axes['An'].plot(XAn, depth_km)
    axes['Jd'].plot(XJd, depth_km)

fig.suptitle('$R_m=${:.1f} km/Myr, $T_m=${:.0f}-{:.0f} °C, $d_g=${}mm'.format(descent_rate/1e3*yr*1e6,T_moho_i-273.15,T_moho_f-273.15,dg/mm,composition.capitalize().replace("_"," ")),y=0.9)
plt.savefig(Path(outputPath,"{}.{}".format("results", "pdf")))
plt.savefig(Path(outputPath,"{}.{}".format("results", "png")))
