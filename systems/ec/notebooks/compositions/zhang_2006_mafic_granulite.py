'''
----------------------------------------

Stable phases at:
                             T(K)     =  1473.00
                             P(bar)   =  20000.0

Phase Compositions (molar  proportions):
                   wt %      vol %     mol %     mol        NA2O     MGO      AL2O3    SIO2     CAO      FEO
 Pl                32.31     38.17     27.13    0.121      0.36874  0.00000  0.63126  2.73747  0.26253  0.00000
 Cpx               24.82     23.39     25.52    0.114      0.13070  0.44817  0.30286  1.82785  0.68491  0.17196
 Gt                34.80     28.69     17.31    0.774E-01  0.00092  1.26321  0.99776  3.00316  0.51475  1.22152
 qtz                8.06      9.75     30.03    0.134      0.00000  0.00000  0.00000  1.00000  0.00000  0.00000

Phase speciation (molar proportions):

 Pl                ab: 0.73747, an: 0.26253
 Cpx               jd: 0.26140, di: 0.34080, hed: 0.17196, cen: 0.05369, cts: 0.17215
 Gt                namj: 0.00092, maj: 0.00132, gr: 0.17158, alm: 0.40717, py: 0.41900

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
    0.2482, # cpx
    0.0000, # opx
    0.0806, # quartz
    0.3231, # plag
    0.3480, # garnet
    0.0000, # kyanite
]

Xik0=[
    [0.34080, 0.17196, 0.05369, 0.17215, 0.26140], # di, hed, cEn, cats, jd
    [1., 0., 0., 0.], # en, fs, mgts, oDi
    [1.], # quartz
    [0.26253, 0.73747], # an, ab
    [0.41900, 0.40717, 0.17158, 0.00132, 0.00092], # py, alm, gr, mgmaj, namaj
    [1.], # kyanite
]