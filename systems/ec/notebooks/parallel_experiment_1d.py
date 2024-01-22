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
reference= 'parallel_experiment_1d'
composition = 'hacker_2003_morb'
rxnName = 'eclogitization_2024_stx21_rx'

# only phases greater than this fraction will be plotted
phasetol = 1.e-2 # 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9
max_steps = 1e3

T_moho_i = 800.
T_moho_f = 800.

shortening_rate_kmMyr = 2.0 # km/Myr = mm/yr
shortening_rate = 2.0*mm/yr # m/yr
L0 = 100. # km

z_moho0 = 30.
z_top0 = 15.

Zs = np.linspace(z_top0,z_moho0,20)

avg_density = 2800 # kg/m3, assumes a linear increase in density from 2600 at the surface to 3000 at the moho

reaction_rate_gcm =  [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-3] # g/cm2/yr
reaction_rate = [r*g/cm/cm/yr for r in reaction_rate_gcm]# kg/m2/yr

# for now, just use 1 rate
reaction_rate = reaction_rate[-2]

print("rxn rate = {} kg/m2/yr".format(reaction_rate))

show_equilibrium = False

# number of processes
processes = mp.cpu_count()

# ------------------------------------------

# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("composition")
    args = parser.parse_args()

    if args.composition is not None:
        composition = args.composition

#====================================================

# Damkhoeler
Da = reaction_rate / (avg_density * shortening_rate) # ????
print("Da={}".format(Da))
# end time of reactions
dt = 100e3 # yrs
total_t = 1e6 # yrs
n_steps = total_t/dt

ts = np.linspace(0, total_t, int(n_steps))

mod = importlib.import_module("compositions."+composition)
Cik0, Xik0, mi0, phii0, phase_names, endmember_names = [getattr(mod,a,None) for a in ['Cik0', 'Xik0', 'mi0','phii0', 'phase_names', 'endmember_names']]

rxn = EcModel.get_reaction(rxnName)

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

rho_final = np.zeros(Zs.shape+(int(n_steps),))
z_final = np.zeros(rho_final.shape)
t_final = np.zeros(rho_final.shape)
# function to run in parallel
# runs for each individual point in the profile
def run_experiment(args):

    zi, z_km = args
    z0 = z_km
    P = 0.028*z_km # GPa
    T = 800+273. # assume constant for now
    _Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
    _mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0

    # Equilibrate model at initial (T0, P0)

    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T,GPa2Bar(P),_mi0,_Cik0,1e5,Da=1e4,eps=eps)

    # get "equilibrated" composition
    mi1 = ode.sol.y[:ode.I, -1]
    Cik1 = ode.sol.y[ode.I:ode.I+ode.K, -1]
    z_top = z_top0
    z_moho = z_moho0
    rhot = np.zeros(ts.shape)
    zt = np.zeros(ts.shape)

    for n,t in enumerate(ts):
        z = z0 + z0/L0*shortening_rate*t
        z_moho = z_moho0 + z_moho0/L0*shortening_rate*t
        z_top = z_top0 + z_top0/L0*shortening_rate*t
        P = 0.028*z_km
        T = 400. + (z-z_top)/(z_moho - z_top)*400 + 273
        ode = ScipyPDReactiveODE(rxn)
        ode.solve(T,GPa2Bar(P),mi1,Cik1,dt,Da=Da,eps=eps,method="BDF_mcm",max_steps=max_steps)
        
        Cik1 = ode.sol.y[ode.I:ode.I+ode.K,-1]
        mi1 = ode.sol.y[:ode.I,-1] # -1 = final time step

        rhot[n] = ode.final_rho()
        zt[n] = z

    return rhot, zt, zi

# run a time series for each depth point
with Pool(processes, maxtasksperchild=12) as pool:
    sols = pool.map(run_experiment, enumerate(Zs))

# loop over depths
for sol in sols:
    rhos, zs, zi = sol
    rho_final[zi,:] = rhos
    z_final[zi,:] = zs
    t_final[zi,:] = ts

outputPath = Path("figs",reference,composition,rxnName)
outputPath.mkdir(parents=True, exist_ok=True)

fig = plt.figure(figsize=(12,14))
ax = fig.add_subplot(1,1,1)
plt.ylim([0,100])
ax.invert_yaxis()

plt.scatter(t_final, z_final, c=rho_final, s=100,alpha=0.75)
plt.show()