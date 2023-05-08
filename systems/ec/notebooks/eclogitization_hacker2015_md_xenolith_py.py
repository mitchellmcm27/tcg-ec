from mcm import EcModel
import matplotlib.pyplot as plt
from pathlib import Path

reference= 'eclogitization_hacker2015_md_xenolith'
rxnName = 'eclogitization_agu5_slb_rx'

outputPath = Path("figs",rxnName,reference)
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



model = EcModel(
    reference,
    rxnName,
    mi0=[
        0.1470, # cpx
        0.2920, # opx
        0.0500, # quartz
        0.5560, # plag
        0.0, # garnet
        0.0, # kyanite
    ],
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
    ],
    P0=1.3,
    T0=1173,
    nP=30,
    nT=30
)

rxn, grid, ode, bdfdiag = model.run(reload=False,save=True,plot_phases=False,end_t=1e4)
plt.savefig(Path(outputPath,"density.png"))
bdfdiag.plot_phases()
plt.savefig(Path(outputPath,"phases.png"))
