from mcm import EcModel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool
import multiprocessing as mp

reference= 'parallel_profile_hacker2015_md_xenolith'

### ------------ INPUTS -------------------

# number of x-nodes
nT = 100

# end time of reactions
end_t = 1e6

# which reaction to use
rxnName = 'eclogitization_agu5_stx21_rx'

# only phases greater than this fraction will be plotted
phasetol = 1.e-3 # 1.e-2

# Damkhoeler number
Da = 1.0 # 1.0
# regularization parameter for compositions
eps = 1.e-3 # 1.e-2

Pmin, Pmax = [2.5, 0.5]
Tmin, Tmax = [773., 1273.]

# number of processes
processes = 20 # mp.cpu_count()

# ------------------------------------------

# Perple_X output...
''' Stixrude 2021

Stable phases at:
                             T(K)     =  1273.00
                             P(bar)   =  5000.00

Phase Compositions (molar  proportions):
                   wt %      vol %     mol %     mol        NA2O     MGO      AL2O3    SIO2     CAO      FEO
 Pl                55.60     61.47     50.34    0.201      0.20011  0.00000  0.79989  2.40023  0.59977  0.00000
 Cpx               14.70     13.10     16.10    0.643E-01  0.02656  0.67820  0.03646  1.99010  0.90294  0.30271
 Opx               29.20     24.85     31.52    0.126      0.00000  1.09205  0.02807  1.97193  0.03959  0.84029
 qtz                0.50      0.58      2.04    0.816E-02  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000

Phase speciation (molar proportions):

 Pl                ab: 0.40023, an: 0.59977
 Cpx               jd: 0.05313, di: 0.59033, hed: 0.30271, cen: 0.04393, cts: 0.00990
 Opx               odi: 0.03959, en: 0.51220, fs: 0.42014, ts: 0.02807

 
'''


phase_list = [
    'Clinopyroxene',
    'Orthopyroxene',
    'Quartz',
    'Feldspar', 
    'Garnet', 
    'Kyanite',
]

em_list = [
    'Diopside', 'Hedenbergite', 'Clinoenstatite', 'CaTschermaks', 'Jadeite',
    'Enstatite', 'Ferrosilite', 'MgTschermaks', 'OrthoDiopside',
    'Quartz',
    'Anorthite','Albite',
    'Pyrope', 'Almandine', 'Grossular', 'MgMajorite', 'NaMajorite',
    'Kyanite'
]

# mass fractions of the phases
## Grt-Opx-Cpx granulite
## given as volume fractions
phii0 = [
    0.1310, # cpx
    0.2485, # opx
    0.00580, # quartz
    0.6147, # plag
    0.0, # garnet
    0.0, # kyanite
 ]

mi0 = [
    0.1470,
    0.2920,
    0.0050,
    0.5560,
    0.0,
    0.0
]

Xik0 = [
    [0.59033, 0.30271, 0.04393, 0.00990, 0.05313], # di, hed, *cEn, *cats, jd
    [0.51220, 0.42014, 0.02807, 0.03959], # en, fs, *mgts, *oDi
    [1.], # quartz
    [0.59977, 0.40023], # an, ab
    [0.39681, 0.42983, 0.17322, 0.0000, 0.0000], # py, alm, gr, *mgmaj, *namaj
    [1.], # kyanite
]

# move cEn to oEn
Xik0[1][1] += Xik0[0][2]
Xik0[0][2] = 0.0
# move oDi to di
Xik0[0][0] += Xik0[1][3]
Xik0[1][3] = 0.0

# regularize 3-component garnet
g3 = (1-(Xik0[4][0]+Xik0[4][1]+Xik0[4][2]))/3.0
Xik0[4][0] += g3
Xik0[4][1] += g3
Xik0[4][2] += g3
Xik0[4][3] = 0.0
Xik0[4][4] = 0.0

rxn = EcModel.get_reaction(rxnName)

def x2c(rxn, Xik0):
    return np.asarray([c for (i, ph) in enumerate(rxn.phases()) for c in ph.x_to_c(Xik0[i])])

Cik0 = x2c(rxn, Xik0)

T_range = np.linspace(Tmin, Tmax, nT)
P_range = np.linspace(Pmin, Pmax, nT)

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
    end = 1e2 if(P>2.25) and (T<800) else end_t
    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T,GPa2Bar(P),mi0,Cik0,end,Da=Da,eps=eps)
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
    print("({},{})".format(T,P))
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
plt.show(block=False)

fig = plt.figure(figsize=(12,12))
axi = fig.add_subplot(1,1,1)
for i, phase in enumerate(rxn.phases()):    
    plt.plot(T_range,phi_final[:,i],':' if i>9 else '-')

plt.ylim([0.005, 0.615])
plt.xlim([Tmin,Tmax])
plt.xticks([873, 973, 1073, 1173, 1273])
plt.yticks([0.005,0.123,0.246,0.369,0.492,0.615])
plt.legend(phase_list)
plt.show(block=False)

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
    "Anorthite": "--",
    "Albite": "--",
    'Pyrope': "--", 
    'Almandine':"--", 
    'Grossular':"--", 
    'MgMajorite':"--", 
    'NaMajorite':"--",
    }

for k in range(K):
    name = em_list[k]
    line_style = line_style_by_endmember[name]
    plt.plot(T_range,Cik_final[:,k], line_style)

plt.ylim([0,1])
plt.xlim([Tmin,Tmax])
plt.legend(em_list)
plt.show(block=False)

plt.show()




    

