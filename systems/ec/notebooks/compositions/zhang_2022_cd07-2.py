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

Xik0=[
    [0.33626, 0.24509, 0.01950, 0.04826, 0.35090], # di, hed, *cEn, *cats, jd
    [1., 0., 0., 0.], # en, fs, *mgts, *oDi
    [1.], # quartz
    [0.26253, 0.73747], # an, ab
    [0.28309, 0.57640, 0.13924, 0.00064, 0.00063], # py, alm, gr, *mgmaj, *namaj
    [1.], # kyanite
]