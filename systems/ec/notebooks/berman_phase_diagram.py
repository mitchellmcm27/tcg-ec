import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
#from __future__ import annotations # Enable Python 4 type hints in Python 3
from thermoengine.equilibrate import PhaseLibrary, GibbsMinimizer, System
import thermoengine as thermo
from thermoengine.const import units
from thermoengine.core import UnorderedList
import os
import contextlib
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import math
import matplotlib.colors as mcolors

from boustrophedon import boustrophedon

def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)
    return wrapper

@supress_stdout
def update_system(sys,T,P, comp=None, phase_library=None):
    if sys is None:
        assert(comp is not None)
        assert(phase_library is not None)
        sys = System(
            T=800,
            P=0.25*units.GPA,
            comp=comp,
            options={'grid_spacing':0.2},
            phase_library=phase_library,
        )
    sys.update(T=T, P=P)
    return sys

def run(phase_symbols = None, oxides = None, T=None, P=None, name=None):

    phase_symbols_berman = phase_symbols if phase_symbols is not None else [
        "Qz", "Coe",
        "Grt",  
        'Di','Jd','cEn', #'Cpx'
        'En', 'Fs',
        "Ky",
        "Ab", "An",
        "Rt",
        "Mc"
    ]

    assert(name is not None)
    assert(P is not None)
    assert(T is not None)
    assert(oxides is not None)

    
    Tmin = np.amin(T)
    Tmax = np.amax(T)
    Pmin = np.amin(P)
    Pmax = np.amax(P)

    oxide_comp = thermo.OxideWtComp(**oxides)

    db = thermo.model.Database(database='Berman')
    phases = db.get_phases(phase_symbols_berman)
    phase_library = PhaseLibrary(phases)

    system = update_system(
            None,
            T=800,
            P=0.25*units.GPA,
            comp=oxide_comp,
            phase_library=phase_library,
        )

    Ts, Ps = np.meshgrid(T,P)

    comps = np.empty(Ts.shape, dtype='U32')
    energy = np.empty(Ts.shape, dtype='float')
    assemblages = np.empty(Ts.shape, dtype='object')
    
    for j,t in enumerate(T):
        for i,p in enumerate(P):
            update_system(system, t,p)
            names = system.stable_assemblage.names[system.stable_assemblage.amounts>0.001]
            hash = " ".join(sorted(names))
            comps[i][j] = hash
            energy[i][j] = system.stable_assemblage.total_energy
            assemblages[i][j] = system.stable_assemblage

    try:
        plot_result(Ts,Ps,assemblages,name)
    except Exception as e: print(e)

    return [Ts, Ps, assemblages]

def plot_result(Ts,Ps,assemblages,name):
    
    get_hash = lambda assemblage : " ".join(sorted(assemblage.names[assemblage.amounts>0.001]))
    comps = np.vectorize(get_hash)(assemblages)

    Tmin = np.amin(Ts)
    Tmax = np.amax(Ts)
    Pmin = np.amin(Ps)
    Pmax = np.amax(Ps)

    categories = np.unique(comps).flatten()
    print(categories)

    colormap = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20)), plt.cm.tab20c(np.linspace(0, 1, 20))))

    colorLookup = dict([(key, value) for key, value in zip(categories, colormap)])
    print(colorLookup.keys())

    colors = [colorLookup[c] for c in comps.flatten()]

    shape = (comps.shape[0], comps.shape[1], 4)
    colorgrid = np.reshape(colors, shape)

    fig = plt.figure(figsize=(10,8))
    ax1 = plt.subplot(4,8,(1,7))
    for q,cat in enumerate(categories):
        indx = comps.flatten() == cat
        ax1.scatter(Ts.flatten()[indx], Ps.flatten()[indx], label=cat, color=colorLookup[cat])
    plt.xlim([Tmin,Tmax])
    plt.ylim([Pmin,Pmax])
    plt.ylabel("Pressure (Pa)")

    ax1.legend(loc='center left')

    ax2 = plt.subplot(4,8,(1,15))
    ax2.imshow(colorgrid, origin='lower', extent=(Tmin,Tmax,Pmin/units.GPA,Pmax/units.GPA), aspect="auto")
    plt.ylabel("Pressure (GPa)")
    plt.xlabel("Temperature (K)")

    ax3 = plt.subplot(4,8,(4,16))
    ax3.legend(*ax1.get_legend_handles_labels(), loc='center left')
    ax3.axis('off')

    ax4 = plt.subplot(4,8,(1,15))
    im=ax2.imshow(n_phases, origin='lower', extent=(Tmin,Tmax,Pmin/units.GPA,Pmax/units.GPA), aspect="auto", cmap="Blues")
    plt.ylabel("Pressure (GPa)")
    plt.xlabel("Temperature (K)")
    plt.colorbar(im)

    fig.suptitle(name)
    plt.savefig(name + "_"+ dt.datetime.now().isoformat() + ".png",facecolor='white', transparent=False)
            
