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
from scipy import ndimage
from scipy.interpolate import griddata
import from_perplex as pp

### ======================= INPUTS ============================
reference = "parallel_pd"
composition = "hacker_2003_morb"
rxnName = "eclogitization_agu10_stx21_rx"

# number of nodes in each dimension (T,P) = (x,y)
nP = 60
nT = 60

# end time of reactions
end_t = 1e5
# limit the maximum number of steps taken by the solver
# defaults to infinity, which can exhaust memory if solver doesn't converge
max_steps = 4e3

# TODO: these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9

Da = 1.0 # Damkhoeler num, default 1.0

# regularization parameter for compositions
# sets the minimum value for Cik
eps = 1.e-5 # default, 1.e-2

Pmin, Pmax = [0.5, 2.5] # Gpa
Tmin, Tmax = [773., 1273.] # K

# number of processes
processes = mp.cpu_count()

# Plotting options
phasetol = 1.e-2 # # phases less than this won't plot, default 1.e-2
density_levels = np.arange(28.5, 37.6, 0.1)
density_ticks = np.asarray([28.,30.,32.,34.,36.,38.])
highlight_densities = np.asarray([33.30] )# /10
density_cmap = plt.get_cmap("Spectral_r")
diff_cmap = plt.get_cmap("bwr_r")
stime_cmap = plt.get_cmap("bwr")
variance_cmap = plt.get_cmap("Blues")
contour_kwargs = {
    "colors": "black",
    "linewidths": 1,
    "linestyles": "solid",
    "alpha": 0.25
}
imshow_kwargs = {
    "extent": (Tmin,Tmax,Pmin,Pmax),
    "origin":"lower",
    "aspect":"auto",
}

fig_filetypes = ["png", "pdf"]

# ===================================================


# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("composition")
    parser.add_argument("-l", "--lowres", action="store_true")
    args = parser.parse_args()

    if args.composition is not None:
        composition = args.composition
    if args.lowres:
        nP = int(nP/2)
        nT = int(nT/2)
        end_t = 1e2

#====================================================

# get reaction object
rxn = EcModel.get_reaction(rxnName)

# import initial compositions from composition file
module = importlib.import_module("compositions."+composition)
Cik0, Xik0, mi0, phii0, phase_names, endmember_names = [
    getattr(module,a,None) for a in [
        "Cik0", 
        "Xik0", 
        "mi0",
        "phii0", 
        "phase_names", 
        "endmember_names"
        ]
    ]

def x2c(rxn, Xik0):
    return np.asarray([c for (i, ph) in enumerate(rxn.phases()) for c in ph.x_to_c(Xik0[i])])
def phi2m(rxn, phii0, Cik0, T=900.,p=10000.):
    densities = []
    C = rxn.zero_C()
    Ki = 0
    for i,ph in enumerate(rxn.phases()):
        n = len(ph.endmembers())
        C[i] = Cik0[Ki:Ki+n]
        Ki = Ki+n
    # regularize C
    C = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    C = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    densities = [ph.rho(T, p, C[i]) for i,ph in enumerate(rxn.phases())]
    mass_tot = np.sum(np.asarray(densities) * np.asarray(phii0))
    mi0 = np.asarray([v*densities[i]/mass_tot for (i, v) in enumerate(phii0)])
    return mi0

Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
mi0 = phi2m(rxn, phii0, Cik0) if mi0 is None else mi0

T_range = np.linspace(Tmin, Tmax, nT)
P_range = np.linspace(Pmin, Pmax, nP)
P_g, T_g  = np.meshgrid(P_range, T_range, indexing="ij")
rho_g = np.zeros(P_g.shape)
stime_g = np.zeros(P_g.shape)
variance_g = np.zeros(P_g.shape)
phases_g = [["" for j in range(nT)] for i in range(nP)]

# function to run in parallel
# solves an ODE for each point in the grid given i, j indices
def task(args):
    i, j = args
    P = P_range[i]
    T = T_range[j]

    ode = ScipyPDReactiveODE(rxn)

    ode.solve(T,GPa2Bar(P),mi0,Cik0,end_t,Da=Da,eps=eps,rtol=rtol,atol=atol,method="BDF_mcm",max_steps=max_steps)
    odephasenames, phaseabbrev = ode.final_phases(phasetol)
    phases = "+".join(phaseabbrev) if ode.sol.status == 0 else ""
    rho = ode.final_rho() if ode.sol.status == 0 else np.nan
    stime = ode.stime
    variance = odephasenames.size
    return rho, phases, stime, variance, i, j

# Solve individual points in parallel
with Pool(processes,maxtasksperchild=12) as pool:
    arglist = [[i,j] for i,P in enumerate(P_range) for j,T in enumerate(T_range)]
    sols = pool.map(task, arglist)

# Collection individual points back into grid
for rho, phases, stime, variance, i , j in sols:
    rho_g[i][j] = rho
    phases_g[i][j] = phases
    stime_g[i][j] = stime
    variance_g[i][j] = variance

# Interplate NaNs in density
x,y = np.indices(rho_g.shape)
rhonan = np.isnan(rho_g)
rho_g[np.isnan(rho_g)] = griddata((x[~rhonan], y[~rhonan]),rho_g[~rhonan],(x[rhonan], y[rhonan]))

# Calculate unique phase assemblage indicies
uniquestrs = sorted(list(set([phstr for phstrl in phases_g for phstr in phstrl if phstr != ""])))
def index(ls,v):
    try:
        i = ls.index(v)
    except ValueError:
        i = -1
    return i
phaseis_g = np.asarray([[index(uniquestrs,phstr) for phstr in phasestrr]  for phasestrr in phases_g])
# Get pyrolite density grid
interp = pp.get_rho_interpolator("data/xu_2008_pyrolite.tab")
rho_pyrolite_g = interp(T_g, GPa2Bar(P_g))

# Create folders if needed
outputPath = Path("figs",reference,composition,rxnName)
outputPath.mkdir(parents=True, exist_ok=True)

def save_current_fig_as(name):
    _name = name[1:] if name[0]=="." else name
    for ext in fig_filetypes:
        _ext = ext[1:] if ext[0]=="." else ext
        plt.savefig(Path(outputPath,"{}.{}".format(_name, _ext)))

def plot_phase_labels(ax):
    fxgrid = T_g.flatten()
    fygrid = P_g.flatten()
    fphaseis = phaseis_g.flatten()

    for i,ph in enumerate(uniquestrs):
        if i in fphaseis:
            indices = fphaseis==i
            ax.text(np.mean(fxgrid[indices]), np.mean(fygrid[indices]),ph,fontsize=12,ha='center',va='center',weight="bold")

# Plot density

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)

s = axi.contourf(T_g,P_g,rho_g,levels=density_levels,alpha=0.75,cmap=density_cmap)
axi.contour(T_g,P_g,rho_g,levels=density_levels,alpha=1,cmap=density_cmap)
axi.contour(T_g,P_g,rho_g,levels=highlight_densities,alpha=1,colors="black")
plt.xlim([Tmin, Tmax])
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location="left")
plt.suptitle(composition)
save_current_fig_as("density")

# Plot phases

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
scs = []
for i,ph in enumerate(uniquestrs):
    scs.append(axi.scatter(T_g[phaseis_g==i],P_g[phaseis_g==i],alpha=0.35,label=ph,s=50))
fig.legend(handles=scs,loc="center left")
plot_phase_labels(axi)
plt.suptitle(composition)


# Plot variance

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
max_variance = np.nanmax(variance_g)
min_variance = np.nanmin(variance_g)
s = axi.scatter(T_g,P_g,c=variance_g,s=50,cmap=variance_cmap,alpha=0.5,vmin=min_variance-1,vmax=max_variance+1)
plot_phase_labels(axi)
plt.suptitle(composition)
save_current_fig_as("phase-variance")

# Plot solution time

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
#s = plt.imshow(stime_g,cmap=stime_cmap,norm=mpl.colors.LogNorm(), **imshow_kwargs)
s = plt.scatter(T_g,P_g,c=stime_g,s=100,alpha=0.75,cmap=stime_cmap,norm=mpl.colors.LogNorm())
fig.colorbar(s,location="left",ax=axi)
plt.suptitle(composition)
save_current_fig_as("stime")

# Plot comparison with Perple_X density
interp = pp.get_rho_interpolator("data/"+composition+".tab")
if (interp is not None):

    interp = pp.get_rho_interpolator("data/hacker_2015_md_xenolith.tab")
    rho_eq_g = interp(T_g, GPa2Bar(P_g))

    fig = plt.figure(figsize=(12,6))
    
    axi = fig.add_subplot(1,3,1)

    s = axi.contourf(T_g, P_g, rho_g, levels=density_levels, alpha=0.75, cmap=density_cmap)
    axi.contour(T_g, P_g, rho_g, levels=density_levels, alpha=1, cmap=density_cmap)
    plt.xlim([Tmin, Tmax])
    plt.colorbar(mappable=s, location="bottom", ticks=density_ticks, label="10$^2$ kg/m$^3$")
    plt.gca().set_title("Reactive density (TCG)")

    axi = fig.add_subplot(1,3,2)
    s = axi.contourf(T_g, P_g, rho_eq_g, levels=density_levels, alpha=0.75, cmap=density_cmap)
    axi.contour(T_g, P_g, rho_eq_g, levels=density_levels, alpha=1, cmap=density_cmap)
    plt.xlim([Tmin, Tmax])
    plt.ylim([Pmin,Pmax])
    plt.colorbar(mappable=s, location="bottom",ticks=density_ticks, label="10$^2$ kg/m$^3$")
    plt.gca().set_title("Equilibrium density (Perple_X)")

    axi = fig.add_subplot(1,3,3)
    levels = np.arange(-2, 2.2, 0.2)
    diff = (rho_g-rho_eq_g)
    s=axi.imshow(diff,cmap=diff_cmap,vmin=-2,vmax=2, **imshow_kwargs)
    axi.contour(T_g ,P_g, diff, levels=levels,**contour_kwargs)
    plt.colorbar(mappable=s, location="bottom",ticks=[-2.,-1.,0.,1.,2.], label="10$^2$ kg/m$^3$")
    plt.gca().set_title("Diff. (TCG minus PX)")
    plt.suptitle(composition)
    save_current_fig_as("density_comparison")

# Plot comparison with pyrolite

fig = plt.figure(figsize=(12,6))

axi = fig.add_subplot(1,3,1)
s = axi.contourf(T_g, P_g, rho_g, levels=density_levels, alpha=0.75, cmap=density_cmap)
axi.contour(T_g, P_g, rho_g, levels=density_levels, alpha=1, cmap=density_cmap)
plt.xlim([Tmin, Tmax])
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location="bottom", ticks=density_ticks, label="10$^2$ kg/m$^3$")
plt.gca().set_title("Reactive density (TCG)")

axi = fig.add_subplot(1,3,2)
s = axi.contourf(T_g, P_g, rho_pyrolite_g, levels=density_levels, alpha=0.75, cmap=density_cmap)
axi.contour(T_g, P_g, rho_pyrolite_g, levels=density_levels, alpha=1, cmap=density_cmap)
plt.xlim([Tmin, Tmax])
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location="bottom", ticks=density_ticks, label="10$^2$ kg/m$^3$")
plt.gca().set_title("Pyrolite density (Xu et al. 2008)")

axi = fig.add_subplot(1,3,3)
diff = (rho_g-rho_pyrolite_g)*100 # g/cm3
maxval = np.ceil(np.nanmax(np.absolute(diff))/100)*100+50
levels = np.arange(-maxval,maxval+50,50)
s=axi.imshow(diff,cmap=diff_cmap,vmin=-maxval,vmax=maxval, **imshow_kwargs)
axi.contour(T_g ,P_g, diff, levels=levels,**contour_kwargs)
plt.colorbar(mappable=s, location="bottom", label="g/cm$^3$")
plt.gca().set_title("Density above pyrolite")
plt.suptitle(composition)
save_current_fig_as("density_contrast")
