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

'''
----------------------------------------

Stable phases at:
                             T(K)     =  1073.15    
                             P(bar)   =  10000.0    

Phase Compositions (weight percentages):
                   wt %      vol %     mol %     mol        NA2O     MGO      AL2O3    SIO2     CAO      FEO  
 Pl                46.68     52.72     41.71    0.172        5.094    0.000   29.232   54.206   11.468    0.000
 Cpx               18.84     17.32     20.68    0.854E-01    1.336   13.379    2.704   54.161   22.398    6.022
 Opx               20.30     17.87     22.09    0.912E-01    0.000   22.979    1.694   53.001    0.308   22.018
 Gt                11.92      9.43      6.42    0.265E-01    0.001   12.000   22.649   40.057    4.307   20.986
 qtz                2.26      2.66      9.10    0.376E-01    0.000    0.000    0.000  100.000    0.000    0.000

Phase speciation (molar proportions):

 Pl                ab: 0.44562, an: 0.55438
 Cpx               jd: 0.09512, di: 0.68544, hed: 0.18496, cen: 0.02353, cts: 0.01096
 Opx               odi: 0.01223, en: 0.60978, fs: 0.34101, ts: 0.03697
 Gt                namj: 0.00005, maj: 0.00025, gr: 0.11521, alm: 0.43819, py: 0.44630
 '''
# mass fractions of the phases

## Eclogite facies
mi0 = [
    0.1884, # cpx
    0.2030, # opx
    0.0226, # quartz
    0.4668, # feldspar
    0.1192, # garnet
    0.000, # kyanite
]

Xik0 = [
    [0.68544, 0.18496, 0.02353, 0.01096, 0.09512], # di, hed, cEn, cats, jd
    [0.60978, 0.34101, 0.03697, 0.01223], # en, fs, mats, oDi
    [1.], # quartz
    [0.55438, 0.44562], # an, ab
    [0.44630, 0.43819, 0.11521, 0.00025, 0.00005], # py, alm, gr, mgmaj, namaj
    [1.], # ky
]

# regularize 3-component garnet
g3 = (1-(Xik0[4][0]+Xik0[4][1]+Xik0[4][2]))/3.0
Xik0[4][0] += g3
Xik0[4][1] += g3
Xik0[4][2] += g3
Xik0[4][3] = 0.0
Xik0[4][4] = 0.0