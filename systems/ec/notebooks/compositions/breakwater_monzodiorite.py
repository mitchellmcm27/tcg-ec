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
    0.16, # cpx
    0.14, # opx
    0.0, # quartz
    0.57+0.05, # plag + k-spar
    0.0, # garnet
    0.08, # kyanite
]

xCpx =  0.229 # Fe/(Fe+Mg)
jCpx =  0.179 # Na/(Ca+Na)
xGrt =  0.588 # Fe/(Fe+Mg)
zGrt =  0.157 # Ca/(Ca+Mg+Fe[2+])
xOpx =  0.364 # Fe/(Fe+Mg)
yOpx =  0.090 # Si+Al-2 ??
NaKfs = 0.248 # Na/(Na+Ca+K)
CaKfs = 0.020 # Ca/(Ca+Na+K)
CaPl =  0.391 # Ca/(Ca+Na+K)
KPl =   0.052 # K/(K+Ca+Na)


# assume all feldspar is plag,
# assume plag is binary soln
# convert from ternary feldspar to binary

NaPl = 1 - CaPl - KPl # Na in plag with K
NaPl = NaPl/(NaPl+CaPl) # renormalize Na without K
CaPl = 1 - NaPl # renormalize Ca without K, Na + Ca = 1
Xik0 = [
    #CaMg           CaFe                       NaAl
    [(1-jCpx)*xCpx, (1-jCpx)*(1-xCpx), 0., 0., jCpx], # di, hed, *cEn, *cats, jd

    #Mg      Fe
    [1-xOpx, xOpx, 0., 0.], # en, fs, *mgts, *oDi

    [1.], # quartz

    #Ca,   Na
    [CaPl, NaPl], # an, ab
    
    #Mg3,               Fe3,           Ca3
    [(1-zGrt)*(1-xGrt), (1-zGrt)*xGrt, zGrt, 0., 0.], # py, alm, gr, *mgmaj, *namaj

    [1.], # kyanite
]

