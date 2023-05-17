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
import scipy.ndimage
import from_perplex as pp

### ------------ INPUTS -------------------
reference = 'parallel_pd'
composition = 'hacker_2015_md_xenolith'
rxnName = 'eclogitization_agu5_stx21_rx'

# number of nodes in each dimension
nP = 45
nT = 45

# end time of reactions
end_t = 1e2

# which reaction to use

# only phases greater than this fraction will be plotted
phasetol = 1.e-2 # 1.e-2

# Damkhoeler number
Da = 1.0 # 1.0
# regularization parameter for compositions
eps = 1.e-3 # 1.e-2

Pmin, Pmax = [0.5, 2.5]
Tmin, Tmax = [475., 1373.]

# number of processes
processes = 20 # mp.cpu_count()
# ------------------------------------------

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
    del ode
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
levels = np.arange(25.,38.1,0.1)

# Resample your data grid by a factor of 3 using cubic spline interpolation.
T_zoomed = scipy.ndimage.zoom(T_final, 5)
P_zoomed = scipy.ndimage.zoom(P_final, 5)
rho_zoomed = scipy.ndimage.zoom(rho_final, 5)

s = axi.contourf(T_zoomed,P_zoomed,rho_zoomed,levels=levels,alpha=0.75,cmap=cmap)
axi.contour(T_zoomed,P_zoomed,rho_zoomed,levels=levels,alpha=1,cmap=cmap)
axi.contour(T_zoomed,P_zoomed,rho_zoomed,levels=[33.30],alpha=1,colors='black')
plt.xlim([Tmin, Tmax])
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location='left')
plt.suptitle(composition)
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
plt.suptitle(composition)
plt.savefig(Path(outputPath,'phases.png'))

# Soln time

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
cmap = plt.get_cmap('bwr_r')
s = plt.scatter(T_final,P_final,c=stime_final,s=100,alpha=0.75,cmap=cmap,norm=mpl.colors.LogNorm())
fig.colorbar(s,location='left',ax=axi)
plt.suptitle(composition)
plt.savefig(Path(outputPath,'stime.png'))

if (composition=='hacker_2015_md_xenolith'):

    interp = pp.get_rho_interpolator('data/hacker_2015_md_xenolith.tab')
    rho_eq = interp(T_final, GPa2Bar(P_final))

    fig = plt.figure(figsize=(12,6))
    
    axi = fig.add_subplot(1,3,1)
    cmap = plt.get_cmap('jet')
    levels = np.arange(2.5, 3.81, 0.01)
    s = axi.contourf(T_final, P_final, rho_final/10, levels=levels, alpha=0.75, cmap=cmap)
    axi.contour(T_final, P_final, rho_final/10, levels=levels, alpha=1, cmap=cmap)
    plt.xlim([Tmin, Tmax])
    plt.colorbar(mappable=s, location='bottom', ticks=[2.5,3.0,3.3,3.6])
    plt.gca().set_title('Reactive density (TCG)')

    axi = fig.add_subplot(1,3,2)
    cmap = plt.get_cmap('jet')
    levels = np.arange(2.5,3.81,0.01)
    s = axi.contourf(T_final, P_final, rho_eq/10, levels=levels, alpha=0.75, cmap=cmap)
    axi.contour(T_final, P_final, rho_eq/10, levels=levels, alpha=1, cmap=cmap)
    plt.xlim([Tmin, Tmax])
    plt.ylim([Pmin,Pmax])
    plt.colorbar(mappable=s, location='bottom',ticks=[2.5,3.0,3.3,3.6])
    plt.gca().set_title('Equilibrium density (Perple_X)')

    axi = fig.add_subplot(1,3,3)
    cmap = plt.get_cmap('bwr')
    levels = np.arange(-0.2, 0.22, 0.02)
    diff = (rho_final-rho_eq)/10
    s=axi.imshow(diff,extent=(Tmin,Tmax,Pmin,Pmax),origin="lower",cmap=cmap,aspect="auto",vmin=-.2,vmax=.2)
    axi.contour(T_final ,P_final, diff, alpha=0.1, levels=levels,color='black')
    plt.colorbar(mappable=s, location='bottom',ticks=[-0.2,-0.1,0,0.1,0.2])
    plt.gca().set_title('Diff. (TCG minus PX)')
    plt.suptitle(composition)
    plt.savefig(Path(outputPath,'density_comparison.png'))


interp = pp.get_rho_interpolator('data/xu_2008_pyrolite.tab')
rho_pyrolite = interp(T_final, GPa2Bar(P_final))

fig = plt.figure(figsize=(12,6))

axi = fig.add_subplot(1,3,1)
cmap = plt.get_cmap('jet')
levels = np.arange(2.5, 3.81, 0.01)
s = axi.contourf(T_final, P_final, rho_final/10, levels=levels, alpha=0.75, cmap=cmap)
axi.contour(T_final, P_final, rho_final/10, levels=levels, alpha=1, cmap=cmap)
plt.xlim([Tmin, Tmax])
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location='bottom', ticks=[2.5,3.0,3.3,3.6])
plt.gca().set_title('Reactive density (TCG)')

axi = fig.add_subplot(1,3,2)
cmap = plt.get_cmap('jet')
levels = np.arange(2.5,3.81,0.01)
s = axi.contourf(T_final, P_final, rho_pyrolite/10, levels=levels, alpha=0.75, cmap=cmap)
axi.contour(T_final, P_final, rho_pyrolite/10, levels=levels, alpha=1, cmap=cmap)
plt.xlim([Tmin, Tmax])
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location='bottom', ticks=[2.5,3.0,3.3,3.6])
plt.gca().set_title('Pyrolite density (Xu et al. 2008)')

axi = fig.add_subplot(1,3,3)
cmap = plt.get_cmap('bwr_r')

diff = (rho_final-rho_pyrolite)/10
maxval = np.amax(np.absolute(diff))
levels = np.arange(-maxval,maxval+0.05,0.05)
s=axi.imshow(diff,extent=(Tmin,Tmax,Pmin,Pmax),origin="lower",cmap=cmap,aspect="auto",vmin=-maxval,vmax=maxval)
axi.contour(T_final ,P_final, diff, alpha=0.1, levels=levels,color='black')
plt.colorbar(mappable=s, location='bottom', ticks=[-maxval,0,maxval])
plt.gca().set_title('Density above pyrolite')
plt.suptitle(composition)
plt.savefig(Path(outputPath,'density_contrast.png'))
