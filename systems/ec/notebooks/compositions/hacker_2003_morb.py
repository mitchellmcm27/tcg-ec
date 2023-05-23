
phase_names = [
    'Clinopyroxene',
    'Orthopyroxene',
    'Quartz',
    'Feldspar', 
    'Garnet', 
    'Kyanite',
]

endmember_names = [
    'Diopside', 'Hedenbergite', 'Clinoenstatite', 'CaTschermaks', 'Jadeite',
    'Enstatite', 'Ferrosilite', 'MgTschermaks', 'OrthoDiopside',
    'Quartz',
    'Anorthite','Albite',
    'Pyrope', 'Almandine', 'Grossular', 'MgMajorite', 'NaMajorite',
    'Kyanite',
]

phii0 = [
        0.185, # cpx
        0.185, # opx
        0.00, # quartz
        0.35, # feldspar
        0.28, # garnet
        0.00, # kyanite
    ]


# mass fractions of the end-members
# given in volume percent, assume constant density among EMs?
Cik0= [
        0.28, 0.61, 0., 0., 0.11, # di, hed, *cEn, *cats, jd
        0.5, 0.5, 0., 0., # en, fs, *mats, *oDi
        1., # quartz
        0.43, 0.57, # an, ab
        0.36, 0.46, 0.18, 0., 0., # py, alm, gr, *mgmaj, *namaj
        1., # ky
    ]