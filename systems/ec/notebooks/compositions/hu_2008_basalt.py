'''
----------------------------------------

Stable phases at:
                             T(K)     =  1073.00
                             P(bar)   =  21000.0

Phase Compositions (molar  proportions):
                   wt %      vol %     mol %     mol        NA2O     MGO      AL2O3    SIO2     CAO      FEO
 Cpx               55.44     57.77     57.71     15.8      0.13772  0.59188  0.16165  1.97607  0.71277  0.12053
 Gt                32.32     29.27     16.21     4.44      0.00016  1.25461  0.99965  3.00051  0.58553  1.15972
 ky                 8.38      7.97     11.64     3.19      0.00000  0.00000  1.00000  1.00000  0.00000  0.00000
 qtz                3.86      4.99     14.44     3.96      0.00000  0.00000  0.00000  1.00000  0.00000  0.00000

Phase speciation (molar proportions):

 Cpx               jd: 0.27545, di: 0.56831, hed: 0.12053, cen: 0.01178, cts: 0.02393
 Gt                namj: 0.00016, maj: 0.00019, gr: 0.19518, alm: 0.38657, py: 0.41790
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
    0.5544, # cpx
    0.0000, # opx
    0.0386, # quartz
    0.3019, # plag
    0.3232, # garnet
    0.0838, # kyanite
]


#Pl     ab: 0.40023, an: 0.59977
#Cpx    jd: 0.05313, di: 0.59033, hed: 0.30271, cen: 0.04393, cts: 0.00990
#Opx    odi: 0.03959, en: 0.51220, fs: 0.42014, ts: 0.02807
Xik0=[
    [0.56831, 0.12053, 0.01178, 0.02393, 0.27545], # di, hed, *cEn, *cats, jd
    [1., 0., 0., 0.], # en, fs, *mgts, *oDi
    [1.], # quartz
    [1., 0.], # an, ab
    [0.41790, 0.38657, 0.19518, 0.00019, 0.00016], # py, alm, gr, *mgmaj, *namaj
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