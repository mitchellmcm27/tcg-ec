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
    'Kyanite'
]

mi0 = [
    0.1470, # cpx
    0.2920, # opx
    0.0050, # quartz
    0.5560, # plag
    0.0, # garnet
    0.0, # kyanite
]

Xik0=[
    [0.59033, 0.30271, 0.04393,  0.00990, 0.05313], # di, hed, cEn, cats, jd
    [0.51220, 0.42014, 0.02807, 0.03959], # en, fs, mgts, oDi
    [1.], # quartz
    [0.59977, 0.40023], # an, ab
    [1., 0., 0., 0., 0.], # py, alm, gr, mgmaj, namaj
    [1.], # kyanite
]