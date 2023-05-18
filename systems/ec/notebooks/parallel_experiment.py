from mcm import EcModel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool
import importlib
import from_perplex as pp

yr = 1 # 3.154e7
mm = 1e-3
g = 1e-3
cm = 1e-2

### ------------ INPUTS -------------------
reference= 'parallel_experiment'
composition = 'hacker_2003_morb'
rxnName = 'eclogitization_agu5_stx21_rx'

# which reaction to use
rxnName = 'eclogitization_agu5_stx21_rx'

# only phases greater than this fraction will be plotted
phasetol = 1.e-2 # 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9
max_steps = 1e4

depth_step_size = 10. # meters
moho_i = 30. # km
moho_f = 70. # km
depth_m = np.arange(moho_i*1e3,moho_f*1e3, depth_step_size) # meter steps
depth_km = depth_m/1e3
P_range = 0.027 * depth_km # GPa
T_moho_i = 800.
T_moho_f = 900.
T_range = T_moho_i + (depth_km - depth_km[0])/(depth_km[-1] - depth_km[0])*(T_moho_f - T_moho_i) + 273.# K - increases by 100 deg? 

shortening_rate_kmMyr = 5.0 # km/Myr = mm/yr
L0 = 100. # km
z0 = moho_i # km
descent_rate_kmMyr = z0/L0 * shortening_rate_kmMyr # mm/yr
descent_rate = descent_rate_kmMyr*mm/yr  # = m/yr
avg_density = 2800 # kg/m3, assumes a linear increase in density from 2600 at the surface to 3000 at the moho
mass_tranport_rate = descent_rate * avg_density # kg/m2/yr
print("mass transport rate = {} kg/m2/yr".format(mass_tranport_rate))
reaction_rate_g =  [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3] # g/cm2/yr
reaction_rate = [r * g/cm/cm/yr for r in reaction_rate_g]# kg/m2/yr
print("rxn rate = {} kg/m2/yr".format(reaction_rate))

# number of processes
processes = 20 # mp.cpu_count()

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
Da = [r/mass_tranport_rate for r in reaction_rate]

# end time of reactions
end_t = depth_step_size/descent_rate # seconds

#print("Da = {}".format(Da))
print("dt = {} kyr".format(end_t/1e3))
print("total t = {} Myr".format(end_t * len(depth_km)/1e6))

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

def run_experiment(Da):
    _Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
    _mi0 = phi2m(rxn, phii0, Cik0) if mi0 is None else mi0

    I = len(rxn.phases())
    Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
    K = sum(Kis)

    rho_final = np.zeros(T_range.shape)
    phases_final = ['' for i,_ in enumerate(depth_km)]
    phi_final = np.empty(T_range.shape+(I,))
    Cik_final = np.empty(T_range.shape+(K,))

    # Equilibrate model at initial (T0, P0)

    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T_range[0],GPa2Bar(P_range[0]),_mi0,_Cik0,1e5,Da=1e4,eps=eps)

    # get "equilibrated" composition
    mi1 = ode.sol.y[:ode.I, -1]
    Cik1 = ode.sol.y[ode.I:ode.I+ode.K, -1]

    # function to run in parallel
    for i,P in enumerate(P_range):
        T = T_range[i]
        ode = ScipyPDReactiveODE(rxn)

        ode.solve(T,GPa2Bar(P),mi1,Cik1,end_t,Da=Da,eps=eps,method="BDF_mcm",max_steps=max_steps)
        odephasenames, phaseabbrev = ode.final_phases(phasetol)
        phases = '+'.join(phaseabbrev)
        rho = ode.final_rho()

        Cik = ode.sol.y[ode.I:ode.I+ode.K,-1]
        mi = ode.sol.y[:ode.I,-1] # -1 = final time step
        mi1 = mi
        Cik1 = Cik
        C = ode.reshapeC(Cik)
        Cs = ode.regularizeC(C)
        Cik_reg = [c for arr in Cs for c in arr]
        rhoi = rxn.rho(T, GPa2Bar(P), Cs)
        v = mi/rhoi
        vi  = 1./v.sum()
        phi = vi*mi/rhoi
        rho_final[i]=rho
        phases_final[i]=phases 
        phi_final[i]=phi 
        #mi_final[i]=mi
        Cik_final[i]=Cik_reg

    return rho_final, phases_final, phi_final, Cik_final, Da



# function to run in parallel
_Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
_mi0 = phi2m(rxn, phii0, Cik0) if mi0 is None else mi0
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

    return rho, phases, phi, Cik_reg, i

# run for varying damkhoeler numbers
with Pool(processes) as pool:
    sols = pool.map(run_experiment, Da)

outputPath = Path("figs",reference,composition,rxnName)
outputPath.mkdir(parents=True, exist_ok=True)

# plot 
fig = plt.figure(figsize=(12,14))
plt.gca().invert_yaxis()

n = len(Da)
colormap = plt.cm.get_cmap('turbo')
plt.gca().set_prop_cycle(plt.cycler('color', colormap(np.linspace(0, 1, n))))

for rho, phases, phi, Cik, _Da in sols:
    plt.plot(rho/10, depth_km)

# run for equilibrium
with Pool(processes,maxtasksperchild=12) as pool:
    sols = pool.map(equilibrate, range(len(P_range)))

# plot
rho_equil = np.zeros(T_range.shape)  
for rho, phases, phi, Cik, i in sols:
    rho_equil[i]=rho

plt.plot(rho_equil/10,depth_km,'k--')

ipyrolite = pp.get_rho_interpolator('data/xu_2008_pyrolite.tab')
rho_pyrolite=ipyrolite((T_range, GPa2Bar(P_range)))
plt.plot(rho_pyrolite/10,depth_km, 'r:')

iharzburgite = pp.get_rho_interpolator('data/xu_2008_harzburgite.tab')
rho_harzburgite=iharzburgite((T_range,GPa2Bar(P_range)))
plt.plot(rho_harzburgite/10,depth_km, 'm:')

plt.legend(['$r_0 = ${:.1e} g/cm$^2$/yr'.format(r) for r in reaction_rate_g] + ['equil.'] + ['pyrolite', 'harzburgite'])
plt.ylabel('Depth (km)')
plt.xlabel('Density')
#plt.xlim([2.9,3.6])
fig.suptitle('$R_m=${:.1f} km/Myr, $T_m=${:.0f}-{:.0f} °C, $d_m=${:.0f}-{:.0f} km, {}'.format(descent_rate_kmMyr,T_moho_i, T_moho_f, moho_i,moho_f,composition),y=0.9)

plt.savefig(Path(outputPath,"{}.{}".format("densities", "pdf")))

