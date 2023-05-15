import os, sys
from mcm import EcModel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool
import multiprocessing as mp
import importlib

### ------------ INPUTS -------------------
reference = 'parallel_pd'
composition = 'hacker_2003_morb'
rxnName = 'eclogitization_agu5_stx21_rx'

# number of nodes in each dimension
nP = 30
nT = 30

# end time of reactions
end_t = 1e4

# which reaction to use

# only phases greater than this fraction will be plotted
phasetol = 1.e-2 # 1.e-2

# Damkhoeler number
Da = 1.0 # 1.0
# regularization parameter for compositions
eps = 1.e-3 # 1.e-2

Pmin, Pmax = [0.5, 2.5]
Tmin, Tmax = [773., 1273.]

# number of processes
processes = 20 # mp.cpu_count()
# ------------------------------------------

importPath = Path("compositions", composition)
print(importPath.absolute())
mod = importlib.import_module("compositions."+composition)

Cik0, Xik0, mi0, phii0, phase_names, endmember_names = [getattr(mod,a,None) for a in ['Cik0', 'Xik0', 'mi0','phii0', 'phase_names', 'endmember_names']]

outputPath = Path("figs",reference,composition,rxnName)
outputPath.mkdir(parents=True, exist_ok=True)

T_range = np.linspace(Tmin, Tmax, nT)
P_range = np.linspace(Pmin, Pmax, nP)


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
    
    #print(densities)
    mass = np.sum(np.asarray(densities) * np.asarray(phii0))
    #print(mass)
    
    mi0 = np.asarray([v*densities[i]/mass for (i, v) in enumerate(phii0)])

    return mi0

Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0

mi0 = phi2m(rxn, phii0, Cik0) if mi0 is None else mi0

P_final, T_final  = np.meshgrid(P_range, T_range, indexing='ij')
rho_final = np.zeros(P_final.shape)
stime_final = np.zeros(P_final.shape)
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
    stime = ode.stime
    return rho, phases, stime, i, j

# Map blocks to block-level calculations
with Pool(processes,maxtasksperchild=12) as pool:
    # blocks until all finished
    arglist = [[i,j] for i,P in enumerate(P_range) for j,T in enumerate(T_range)]
    sols = pool.map(task, arglist)

for rho, phases, stime, i , j in sols:
    rho_final[i][j] = rho
    phases_final[i][j] = phases
    stime_final[i][j] = stime

fig = plt.figure(figsize=(12,14))

axi = fig.add_subplot(1,1,1)
cmap = plt.get_cmap('jet')

levels = np.arange(25.,38.,0.1)

plt.scatter(T_final,P_final,c=rho_final,s=100,alpha=0.75,cmap=cmap)
s = axi.contourf(T_final,P_final,rho_final,levels=levels,alpha=0.75,cmap=cmap)
cs = axi.contour(T_final,P_final,rho_final,levels=levels,alpha=0.75,cmap=cmap)

# bold the 3,330 kg/m3 contour
cs.collections[78].set_linewidth(3)        
cs.collections[78].set_color('black')

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

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
scs = []
for i,ph in enumerate(uniquestrs):
    scs.append(axi.scatter(T_final[phaseis_final==i],P_final[phaseis_final==i],alpha=0.5,label=ph,s=50))
fig.legend(handles=scs,loc='center left')
plt.savefig(Path(outputPath,'phases.png'))

# Soln time

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
cmap = plt.get_cmap('bwr')
s = plt.scatter(T_final,P_final,c=stime_final,s=100,alpha=0.75,cmap=cmap,norm=mpl.colors.LogNorm())
fig.colorbar(s,location='left',ax=axi)
plt.savefig(Path(outputPath,'stime.png'))

