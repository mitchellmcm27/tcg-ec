'''
----------------------------------------

Stable phases at:
                             T(K)     =  1173.00
                             P(bar)   =  13000.0

Phase Compositions (molar  proportions):
                   wt %      vol %     mol %     mol        NA2O     MGO      AL2O3    SIO2     CAO      FEO
 Pl                30.19     35.08     28.77    0.107      0.16182  0.00000  0.83818  2.32363  0.67637  0.00000
 Cpx               42.68     41.25     51.28    0.190      0.05942  0.81200  0.09818  1.96124  0.83481  0.07676
 Gt                26.17     22.49     15.84    0.587E-01  0.00010  1.93542  0.99946  3.00064  0.42421  0.64059
 qtz                0.95      1.17      4.12    0.153E-01  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000

Phase speciation (molar proportions):

 Pl                ab: 0.32363, an: 0.67637
 Cpx               jd: 0.11883, di: 0.71929, hed: 0.07676, cen: 0.04635, cts: 0.03876
 Gt                namj: 0.00010, maj: 0.00043, gr: 0.14140, alm: 0.21353, py: 0.64453
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
    0.4268, # cpx
    0.0000, # opx
    0.0095, # quartz
    0.3019, # plag
    0.2617, # garnet
    0.0, # kyanite
]

Xik0=[
    [0.71929, 0.07676, 0.04635, 0.03876, 0.11883], # di, hed, cEn, cats, jd
    [1., 0., 0., 0.], # en, fs, mgts, oDi
    [1.], # quartz
    [0.67637, 0.32363], # an, ab
    [0.64453, 0.21353, 0.14140, 0.00043, 0.00010], # py, alm, gr, *mgmaj, *namaj
    [1.], # kyanite
]