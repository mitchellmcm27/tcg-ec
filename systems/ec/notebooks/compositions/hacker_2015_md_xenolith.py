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


#Pl     ab: 0.40023, an: 0.59977
#Cpx    jd: 0.05313, di: 0.59033, hed: 0.30271, cen: 0.04393, cts: 0.00990
#Opx    odi: 0.03959, en: 0.51220, fs: 0.42014, ts: 0.02807
Xik0=[
    [0.59033, 0.30271, 0.04393,  0.00990, 0.05313], # di, hed, *cEn, *cats, jd
    [0.51220, 0.42014, 0.02807, 0.03959], # en, fs, *mgts, *oDi
    [1.], # quartz
    [0.59977, 0.40023], # an, ab
    [0.39681, 0.42983, 0.17322, 0., 0.], # py, alm, gr, *mgmaj, *namaj
    [1.], # kyanite
]

# move cEn to oEn
Xik0[1][1] += Xik0[0][2]
Xik0[0][2] = 0.0

# move oDi to di
Xik0[0][0] += Xik0[1][3]
Xik0[1][3] = 0.0

# regularize 3-component garnet
g3 = (1-(Xik0[4][0]+Xik0[4][1]+Xik0[4][2]))/3.0
Xik0[4][0] += g3
Xik0[4][1] += g3
Xik0[4][2] += g3
Xik0[4][3] = 0.0
Xik0[4][4] = 0.0