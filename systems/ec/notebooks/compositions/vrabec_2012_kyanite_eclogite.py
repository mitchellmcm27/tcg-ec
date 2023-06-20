'''
----------------------------------------

Stable phases at:
                             T(K)     =  1073.15    
                             P(bar)   =  10000.0    

Phase Compositions (weight percentages):
                   wt %      vol %     mol %     mol        NA2O     MGO      AL2O3    SIO2     CAO      FEO  
 Pl                34.54     39.80     33.52    0.122        3.544    0.000   31.489   50.855   14.112    0.000
 Cpx               40.33     38.64     49.33    0.179        1.515   15.586    3.343   54.926   22.187    2.443
 Opx                2.62      2.48      3.30    0.120E-01    0.000   31.271    2.880   55.411    0.273   10.165
 Gt                22.52     19.07     13.86    0.504E-01    0.001   18.179   23.658   41.839    4.614   11.710

Phase speciation (molar proportions):

 Pl                ab: 0.31246, an: 0.68754
 Cpx               jd: 0.10599, di: 0.76601, hed: 0.07372, cen: 0.03620, cts: 0.01808
 Opx               odi: 0.01023, en: 0.78147, fs: 0.14886, ts: 0.05943
 Gt                namj: 0.00005, maj: 0.00022, gr: 0.11816, alm: 0.23409, py: 0.64748
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
    0.4033, # cpx
    0.0262, # opx
    0.0000, # quartz
    0.3454, # plag
    0.2252, # garnet
    0.0, # kyanite
]

Xik0=[
    [0.76601, 0.07372, 0.03620, 0.01808, 0.10599], # di, hed, cEn, cats, jd
    [0.78147, 0.14886, 0.05943, 0.01023], # en, fs, mgts, oDi
    [1.], # quartz
    [0.68754, 0.31246], # an, ab
    [0.64748, 0.23409, 0.11816, 0.00022, 0.00005], # py, alm, gr, *mgmaj, *namaj
    [1.], # kyanite
]

# regularize 3-component garnet
g3 = (1-(Xik0[4][0]+Xik0[4][1]+Xik0[4][2]))/3.0
Xik0[4][0] += g3
Xik0[4][1] += g3
Xik0[4][2] += g3
Xik0[4][3] = 0.0
Xik0[4][4] = 0.0