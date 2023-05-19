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

phii0 = [
    0.47, # cpx
    0.17, # opx
    0.0, # quartz
    0.22, # feldspar
    0.0, # garnet
    0.10, # kyanite
 ]


xCpx = 0.265 # Fe/(Fe+Mg)
jCpx = 0.115 # Na/(Ca+Na)
xGrt = 0.637 # Fe/(Fe+Mg)
zGrt = 0.187 # Ca/(Ca+Mg+Fe[2+])
xOpx = 0.413 # Fe/(Fe+Mg)
yOpx = 0.099 # Si+Al-2 ??
NaKfs = 0.168 # Na/(Na+Ca+K)
CaKfs = 0.019 #Ca/(Ca+Na+K)
CaPl = 0.571 # Ca/(Ca+Na+K)
KPl = 0.034 # K/(K+Ca+Na)

# convert from ternary feldspar to binary
NaPl = 1 - CaPl - KPl
NaPl = NaPl/(NaPl+CaPl)
CaPl = 1-NaPl

Xik0 = [
    #CaMg           CaFe                       NaAl
    [(1-jCpx)*xCpx, (1-jCpx)*(1-xCpx), 0., 0., jCpx], # di, hed, *cEn, *cats, jd

    #Mg   Fe
    [1-xOpx, xOpx, 0., 0.], # en, fs, *mgts, *oDi

    [1.], # quartz

    #CaAl2, NaAlSi
    [CaPl, NaPl], # an, ab

    #Mg3,  Fe3,  Ca3
    [(1-zGrt)*(1-xGrt), (1-zGrt)*xGrt, zGrt, 0., 0.], # py, alm, gr, *mgmaj, *namaj
    
    [1.], # kyanite
]