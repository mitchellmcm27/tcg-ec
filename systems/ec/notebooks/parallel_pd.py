from mcm import EcModel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool
import multiprocessing as mp

### ------------ INPUTS -------------------
reference = 'parallel_pd_hacker2015_md_xenolith_py_multi'

# number of nodes in each dimension
nP = 30
nT = 30

# end time of reactions
end_t = 1e4

# which reaction to use
rxnName = 'eclogitization_agu5_slb_rx'

# only phases greater than this fraction will be plotted
phasetol = 2.e-2 # 1.e-2

# Damkhoeler number
Da = 1.0 # 1.0
# regularization parameter for compositions
eps = 1.e-2 # 1.e-2

Pmin, Pmax = [0.5, 2.5]
Tmin, Tmax = [773., 1273.]

# number of processes
processes = 20 # mp.cpu_count()

# ------------------------------------------

outputPath = Path("figs",reference,rxnName)
outputPath.mkdir(parents=True, exist_ok=True)

phases = [
    'Clinopyroxene',
    'Orthopyroxene',
    'Quartz',
    'Feldspar', 
    'Garnet', 
    'Kyanite',
]

ems = [
    'Diopside', 'Hedenbergite', 'Clinoenstatite', 'CaTschermaks', 'Jadeite',
    'Enstatite', 'Ferrosilite', 'MgTschermaks', 'OrthoDiopside',
    'Quartz',
    'Anorthite','Albite',
    'Pyrope', 'Almandine', 'Grossular', 'MgMajorite', 'NaMajorite',
    'Kyanite'
]

T_range = np.linspace(Tmin, Tmax, nT)
P_range = np.linspace(Pmin, Pmax, nP)

mi0 = [
    0.1470, # cpx
    0.2920, # opx
    0.0050, # quartz
    0.5560, # plag
    0.0, # garnet
    0.0, # kyanite
]

#Pl     ab: 0.40023, an: 0.59977
#Cpx    jd: 0.05313, di: 0.59033, hed: 0.30271, cen: 0.04393, cts: 0.00990
#Opx    odi: 0.03959, en: 0.51220, fs: 0.42014, ts: 0.02807
Xik0=[
    [0.59033, 0.30271, 0.04393,  0.00990, 0.05313], # di, hed, *cEn, *cats, jd
    [0.51220, 0.42014, 0.02807, 0.03959], # en, fs, *mgts, *oDi
    [1.], # quartz
    [0.59977, 0.40023], # an, ab
    [0.39681, 0.42983, 0.17322, 0., 0.], # py, alm, gr, *mgmaj, *namaj
    [1.], # kyanite
]

print(Xik0[0][2])
# move cEn to oEn
Xik0[1][1] += Xik0[0][2]
Xik0[0][2] = 0.0
print(Xik0[0][2])

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

P_final, T_final  = np.meshgrid(P_range, T_range, indexing='ij')
rho_final = np.zeros(P_final.shape)
phases_final = [['' for j in range(nT)] for i in range(nP)]

# function to run in parallel
def task(args):
    i, j = args
    P = P_range[i]
    T = T_range[j]
    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T,GPa2Bar(P),mi0,Cik0,end_t,Da=Da,eps=eps,rtol=1.e-4,atol=1.e-8)
    odephasenames, phaseabbrev = ode.final_phases(phasetol)
    phases = '+'.join(phaseabbrev)
    rho = ode.final_rho()
    return rho, phases, i, j

# Map blocks to block-level calculations
with Pool(processes,maxtasksperchild=12) as pool:
    # blocks until all finished
    arglist = [[i,j] for i,P in enumerate(P_range) for j,T in enumerate(T_range)]
    sols = pool.map(task, arglist)

for rho, phases, i , j in sols:
    rho_final[i][j] = rho
    phases_final[i][j] = phases

fig = plt.figure(figsize=(12,14))

axi = fig.add_subplot(1,1,1)
cmap = plt.get_cmap('jet')
levels = np.arange(25.,38.,0.1)
plt.scatter(T_final,P_final,c=rho_final,s=100,alpha=0.75,cmap=cmap)
s = axi.contourf(T_final,P_final,rho_final,levels=levels,alpha=0.75,cmap=cmap)
axi.contour(T_final,P_final,rho_final,levels=levels,alpha=0.75,cmap=cmap)
plt.clim([25.,38.])
plt.xlim([Tmin, Tmax])
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location='left')
plt.savefig(Path(outputPath,'density.png'))

# a set of unique strings representing the phase combinations present across the phase diagram
uniquestrs = sorted(list(set([phstr for phstrl in phases_final for phstr in phstrl if phstr != ''])))
def index(ls,v):
    try:
        i = ls.index(v)
    except ValueError:
        i = -1
    return i
# an index into the uniquestr list to give a grid showing which phases are present where
phaseis_final = np.asarray([[index(uniquestrs,phstr) for phstr in phasestrr]  for phasestrr in phases_final])

fig = plt.figure(figsize=(24,14))
axes = fig.subplot_mosaic([['A',[['1','2'],['3','4'],['5','6']]]])
axi = axes['A']
raxes = [axes[repr(i)] for i in range(1,7)]
for axis in raxes: axis.set_visible(False)
scs = []
for i,ph in enumerate(uniquestrs):
    scs.append(axi.scatter(T_final[phaseis_final==i],P_final[phaseis_final==i],alpha=0.5,label=ph,s=50))
fig.legend(handles=scs,loc='center left')
plt.savefig(Path(outputPath,'phases.png'))


