# Note, Fe2O3T was converted to FeO by multipltying by 0.9

'''
----------------------------------------

Stable phases at:
                             T(K)     =  1273.00
                             P(bar)   =  20000.0

Phase Compositions (molar  proportions):
                   wt %      vol %     mol %     mol        NA2O     MGO      AL2O3    SIO2     CAO      FEO
 Cpx               52.35     54.06     50.63    0.230      0.17545  0.37525  0.22371  1.95174  0.62960  0.24509
 Gt                38.65     34.03     17.63    0.800E-01  0.00063  0.85244  0.99873  3.00191  0.41773  1.72919
 qtz                9.01     11.91     31.74    0.144      0.00000  0.00000  0.00000  1.00000  0.00000  0.00000

Phase speciation (molar proportions):

 Cpx               jd: 0.35090, di: 0.33626, hed: 0.24509, cen: 0.01950, cts: 0.04826
 Gt                namj: 0.00063, maj: 0.00064, gr: 0.13924, alm: 0.57640, py: 0.28309
 
 '''

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
    0.5235, # cpx
    0.0000, # opx
    0.0901, # quartz
    0.0000, # plag
    0.3865, # garnet
    0.0000, # kyanite
]


#Pl     ab: 0.40023, an: 0.59977
#Cpx    jd: 0.05313, di: 0.59033, hed: 0.30271, cen: 0.04393, cts: 0.00990
#Opx    odi: 0.03959, en: 0.51220, fs: 0.42014, ts: 0.02807
Xik0=[
    [0.33626, 0.24509, 0.01950, 0.04826, 0.35090], # di, hed, *cEn, *cats, jd
    [1., 0., 0., 0.], # en, fs, *mgts, *oDi
    [1.], # quartz
    [0.26253, 0.73747], # an, ab
    [0.28309, 0.57640, 0.13924, 0.00064, 0.00063], # py, alm, gr, *mgmaj, *namaj
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