from mcm import EcModel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.base import PDReactiveGrid, PDReactiveProfile, PDReactiveGridDiagnostics, PDReactiveProfileDiagnostics
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE

from multiprocessing import Pool

reference= 'eclogitization_hacker2015_md_xenolith_py_multi'

# %% INPUTS

''' multiprocessing
 Break the domain up into blocks for each cpu

   ncols
 +----+----+
 | 1  | 2  |
 |    |    |
 +----+----+ nrows
 | 3  | 4  |
 |    |    |
 +----+----+

'''

# number of y-dir blocks
nrows = 6 
# number of x-dir blocks
ncols = 6 
# total number of processes to use
nblocks = nrows*ncols

n = 5 # number of computational nodes per dimension per block (nxn grid)

# end time of reactions
end_t = 1e4

# which reaction to use
rxnName = 'eclogitization_agu5_slb_rx'
###

# %% Setup

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


Pmin = 0.5
Pmax = 2.5
Tmin = 773.
Tmax = 1273.

rowcol = [[row, col] for row in range(nrows) for col in range(ncols)]

T_range = np.linspace(Tmin, Tmax, ncols*n)
P_range = np.linspace(Pmin, Pmax, nrows*n)

mi0 = [
    0.1470, # cpx
    0.2920, # opx
    0.0500, # quartz
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
    [0.3, 0.3, 0.4, 0., 0.], # py, alm, gr, *mgmaj, *namaj
    [1.], # kyanite
]

rxn = EcModel.get_reaction(rxnName)

def x2c(rxn, Xik0):
    return np.asarray([c for (i, ph) in enumerate(rxn.phases()) for c in ph.x_to_c(Xik0[i])])

Cik0 = x2c(rxn, Xik0)

# %% Solving

# function to run in parallel
def task(rowcol):
    row, col = rowcol

    t0, t1 = [col*n, (col+1)*n]
    p0, p1 = [row*n, (row+1)*n]

    # block-level grid solve
    bdfgrid = PDReactiveGrid()
    bdfgrid.solve(
        rxn, 
        ScipyPDReactiveODE, 
        2, # i0 - doesn't matter because we pass Cik0
        ['T', 'p'], # [x, y] type for calculation
        T_range[t0:t1], 
        P_range[p0:p1],
        end_t,
        Cik0=Cik0, 
        mi0=mi0
    )
    bdfdiag = PDReactiveGridDiagnostics(rxn, bdfgrid)

    # outputs
    rhogrid = bdfdiag.rhogrid()
    if bdfdiag.phaseis is None: bdfdiag.phase_diagnostics()
    phasegrid = bdfdiag.phasestrs.copy()

    print('{}_{}'.format(reference, row*ncols + col))

    del bdfdiag
    del bdfgrid
    return rhogrid, phasegrid

P_final, T_final  = np.meshgrid(P_range, T_range, indexing='ij')
rho_final = np.zeros(P_final.shape)
phases_final = [['' for j in range(ncols*n)] for i in range(nrows*n)]

# Map blocks to block-level calculations
with Pool(maxtasksperchild=1) as pool:
    # blocks until all finished
    sols = pool.map(task, rowcol)

# Solving finished

# %% Plotting
# loop through CPU blocks to copy outputs into big arrays
for i, rc in enumerate(rowcol):
    row, col = rc
    output = sols[i]

    rho, phases = output
    r0, r1 = [n*row, n*(row+1)]
    c0, c1= [n*col, n*(col+1)]

    # copy into numpy arrays
    rho_final[r0:r1,c0:c1] = rho
    # copy into the python list, loop over each point
    for i in range(n):
        for j in range(n):
            phases_final[r0+j][c0+i] = phases[j][i]

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

try:
    plt.colorbar(mappable=s, location='left')
except:
    pass

plt.savefig(Path(outputPath,'density.png'))

# a set of unique strings representing the phase combinations present across the phase diagram
uniquestrs = sorted(list(set([phstr for phstrl in phases_final for phstr in phstrl if phstr != ''])))
print(uniquestrs)

def index(ls,v):
    try:
        i = ls.index(v)
    except ValueError:
        i = -1
    return i

print(phases_final[0][:])
# an index into the uniquestr list to give a grid showing which phases are present where
phaseis_final = np.asarray([[index(uniquestrs,phstr) for phstr in phasestrr]  for phasestrr in phases_final])

fig = plt.figure(figsize=(24,14))
axes = fig.subplot_mosaic([['A',[['1','2'],['3','4'],['5','6']]]])
axi = axes['A']
raxes = [axes[repr(i)] for i in range(1,7)]
for axis in raxes: axis.set_visible(False)

#jgrid, igrid = np.meshgrid(range(len(self.grid.x_range)), range(len(self.grid.y_range)))
scs = []
sis = []
sjs = []
for i,ph in enumerate(uniquestrs):
    sis.append(P_final[phaseis_final==i])
    sjs.append(T_final[phaseis_final==i])
    scs.append(axi.scatter(T_final[phaseis_final==i],P_final[phaseis_final==i],alpha=0.5,label=ph,s=50))
fig.legend(handles=scs,loc='center left')
plt.savefig(Path(outputPath,'phases.png'))


