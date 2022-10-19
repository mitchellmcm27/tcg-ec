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

def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)
    return wrapper

@supress_stdout
def update_system(sys,T,P):
    sys.update(T=T, P=P)
    amts = sys.stable_assemblage.amounts[sys.stable_assemblage.amounts>0.01]
    names = sys.stable_assemblage.names[sys.stable_assemblage.amounts>0.01]
    names_sorted = [x for _, x in sorted(zip(amts, names), reverse=True)]
    amts_sorted = sorted(amts, reverse=True)
    hash = " ".join(names_sorted)
    return hash

def run(phase_symbols = None, oxides = None, T=None, P=None, name=None):

    phase_symbols_berman = phase_symbols if phase_symbols is not None else [
        "Qz", "Coe",
        "Grt",  
        'Di','Jd','cEn', #'Cpx'
        'En', 'Fs', 'En',
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

    system = System(
        T=800,
        P=0.25*units.GPA,
        comp=oxide_comp,
        options={'grid_spacing':0.2},
        phase_library=phase_library
    )

    Ts, Ps = np.meshgrid(T,P)
    comps = np.empty(Ts.shape, dtype='U32')
    for j,t in enumerate(T):
        for i,p in enumerate(P):
            comps[i][j] = update_system(system, t,p)

    categories = np.unique(comps).flatten()
    cmap = plt.get_cmap('tab20')
    colormap = cmap(np.linspace(0, 1, 20))
    colorLookup = dict([(key, value) for key, value in zip(categories, colormap)])

    colors = [colorLookup[c] for c in comps.flatten()]

    shape = (comps.shape[0], comps.shape[1], 4)
    colorgrid = np.reshape(colors, shape)

    fig =plt.figure(figsize=(8,10),facecolor='white')
    ax  =plt.subplot(211)
    for q,cat in enumerate(categories):
        indx = comps.flatten() == cat
        ax.scatter(Ts.flatten()[indx], Ps.flatten()[indx], label=cat, color=colorLookup[cat])
    plt.xlim([Tmin,Tmax])
    plt.ylim([Pmin,Pmax])
    plt.ylabel("Pressure (Pa)")

    leg1 = ax.legend()

    ax2 = plt.subplot(212)
    ax2.imshow(colorgrid, origin='lower', extent=(Tmin,Tmax,Pmin/units.GPA,Pmax/units.GPA), aspect="auto")
    ax2.legend(*ax.get_legend_handles_labels())
    plt.ylabel("Pressure (GPa)")
    plt.xlabel("Temperature (K)")
    plt.title(name)
    plt.savefig(name + "_"+ dt.datetime.now().isoformat() + ".png")

    return [
        Ts, Ps, comps
    ]
            
