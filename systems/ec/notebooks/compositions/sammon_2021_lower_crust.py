phase_names = [
    'Clinopyroxene',
    'Orthopyroxene',
    'Quartz',
    'Feldspar', 
    'Garnet', 
    'Kyanite',
]

endmember_names =[
    'Diopside', 'Hedenbergite', 'Clinoenstatite', 'CaTschermaks', 'Jadeite',
    'Enstatite', 'Ferrosilite', 'MgTschermaks', 'OrthoDiopside',
    'Quartz',
    'Anorthite','Albite',
    'Pyrope', 'Almandine', 'Grossular', 'MgMajorite', 'NaMajorite',
    'Kyanite'
]


# mass fractions of the phases

# initialize with eclogitic compositions
# 40% cpx, 33% garnet
# Mg# should be around 50

## Eclogite facies
mi0_eclogite = [
    0.40, # cpx
    0.00, # opx
    0.20, # quartz
    0.00, # feldspar
    0.33, # garnet
    0.07, # kyanite
]

# Granulite facies
mi0_granulite = [
    0.09, # cpx
    0.24, # opx
    0.08, # quartz
    0.51, # feldspar
    0.0, # garnet
    0.07, # kyanite
]

mi0 = mi0_eclogite

# Cik0 = mass fractions
# Xik0 = mol fractions
# note: * = thermodynamic endmember, set to zero
Xik0 = [
    [0.25, 0.25, 0., 0., 0.5], # di, hed, *cEn, *cats, jd
    [0.5, 0.5, 0., 0.], # en, fs, *mats, *oDi
    [1.], # quartz
    [0.25, 0.75], # an, ab
    [0.4, 0.4, 0.2, 0., 0.], # py, alm, gr, *mgmaj, *namaj
    [1.], # ky
]
