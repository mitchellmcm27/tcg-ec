from mcm import EcModel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
import multiprocessing as mp
from multiprocessing import Pool
import importlib
import from_perplex as pp

### ------------ INPUTS -------------------
reference= 'parallel_profile'
composition = 'hacker_2015_md_xenolith'
rxn_name = 'eclogitization_agu14_stx21_rx'

# number of x-nodes
nT = 100

# end time of reactions
end_t = 1e5

# only phases greater than this fraction will be plotted
phasetol = 1.e-3 # 1.e-2

# Damkhoeler number
Da = 1.0 # 1.0
# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9

# large number
max_steps = 1e5 # 4e3 is reasonable

Pmin, Pmax = [2.5, 0.5]
Tmin, Tmax = [773., 1273.]

# number of processes
processes = mp.cpu_count()

# ------------------------------------------

# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-e", "--end_t")
    parser.add_argument("-r", "--rxn_name")

    args = parser.parse_args()

    if args.composition is not None:
        composition = args.composition
    if args.end_t is not None:
        end_t = args.end_t
    if args.rxn_name is not None:
        rxn_name = args.rxn_name

#====================================================

mod = importlib.import_module("compositions."+composition)
Cik0, Xik0, mi0, phii0, phase_names, endmember_names = [getattr(mod,a,None) for a in ['Cik0', 'Xik0', 'mi0','phii0', 'phase_names', 'endmember_names']]

outputPath = Path("figs",reference,composition,rxn_name)
outputPath.mkdir(parents=True, exist_ok=True)

T_range = np.linspace(Tmin, Tmax, nT)
P_range = np.linspace(Pmin, Pmax, nT)

rxn = EcModel.get_reaction(rxn_name)
print(rxn.report())

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

Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
mi0 = phi2m(rxn, phii0, Cik0) if mi0 is None else mi0

I = len(rxn.phases())
Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(Kis)

rho_final = np.zeros(T_range.shape)
phases_final = ['' for i in range(nT)]
phi_final = np.empty(T_range.shape+(I,))
Cik_final = np.empty(T_range.shape+(K,))

# function to run in parallel
def task(i):
    P = P_range[i]
    T = T_range[i]
    # top-left corner is very slow, decrease end_t
    end = end_t/100. if(P>2.25) and (T<800) else end_t

    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T,GPa2Bar(P),mi0,Cik0,end,Da=Da,eps=eps,method="BDF_mcm",max_steps=max_steps)
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

# Map blocks to block-level calculations
with Pool(processes,maxtasksperchild=12) as pool:
    # blocks until all finished
    sols = pool.map(task, range(len(T_range)))

for rho, phases, phi, mi, Cik, i in sols:
    rho_final[i] = rho
    phases_final[i] = phases
    phi_final[i,:] = phi
    Cik_final[i,:] = Cik

fig = plt.figure(figsize=(12,14))

plt.plot(P_range,rho_final)
plt.savefig(Path(outputPath,'density.png'))

fig = plt.figure(figsize=(12,12))
axi = fig.add_subplot(1,1,1)

df = pp.get_profile_data("data/"+composition+"_profile.tab")

# T(K), P(bar), Pl, Pl, Cpx, Opx, qtz, Gt, ky  
phase_i_to_col_i = {
    3:2,
    "none":3,
    0:4,
    1:5,
    2:6,
    4:7,
    5:8,
}
print(df)

hs = []

for i, phase in enumerate(rxn.phases()):
    h = plt.plot(T_range-273.15,phi_final[:,i],':' if i>9 else '-')
    hs.append(h)
plt.legend(phase_names)
plt.xlim([Tmin-273.15,Tmax-273.15])
plt.xticks([500, 600, 700, 800, 900, 1000])
plt.gca().set_ylim(bottom=0)
if(df is not None):
    for i, phase in enumerate(rxn.phases()):
        col = phase_i_to_col_i[i]
        h = hs[i]
        print(h)
        plt.plot(df["T(K)"]-273.15,df.iloc[:,col]/100, "--", alpha=0.8, linewidth=1, color=h[-1].get_color())

plt.savefig(Path(outputPath,'phases.png'))

fig = plt.figure(figsize=(12,12))
axi = fig.add_subplot(1,1,1)

line_style_by_endmember = {
    "Quartz": "-",
    "Kyanite": "-",
    "Diopside": ":",
    "Hedenbergite": ":",
    "Jadeite": ":",
    "CaTschermaks": ":",
    "Clinoenstatite": ":",
    'Enstatite':"-.", 
    'Ferrosilite':"-.", 
    'MgTschermaks':"-.", 
    'OrthoDiopside':"-.",
    "Anorthite": "-",
    "Albite": "-",
    'Pyrope': "--", 
    'Almandine':"--", 
    'Grossular':"--", 
    'MgMajorite':"--", 
    'NaMajorite':"--",
    }

for k in range(K):
    name = endmember_names[k]
    line_style = line_style_by_endmember[name]
    plt.plot(T_range-273.15,Cik_final[:,k], line_style)

plt.ylim([0,1])
plt.xlim([Tmin-273.15,Tmax-273.15])
plt.legend(endmember_names)

plt.savefig(Path(outputPath,'endmembers.png'))