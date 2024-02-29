# what about the following SLB phases:
# spinel
# coesite
# stishovite
# wuestite--periclase SS (magnesiowuestite)? (FeO-MgO)

# it would also be easy(?) to add pure-endmember phases from the HGP database

an = ['Feldspar', 'Anorthite']
ab = ['Feldspar', 'Albite']

fs = ['Orthopyroxene', 'Ferrosilite']
en = ['Orthopyroxene', 'Enstatite']
mgts = ['Orthopyroxene', 'MgTschermaks']
odi = ['Orthopyroxene', 'OrthoDiopside']

hed = ['Clinopyroxene', 'Hedenbergite']
di = ['Clinopyroxene', 'Diopside']
jd = ['Clinopyroxene', 'Jadeite']
cats = ['Clinopyroxene', 'CaTschermaks']
cen = ['Clinopyroxene', 'Clinoenstatite']

alm = ['Garnet', 'Almandine']
pyp = ['Garnet', 'Pyrope']
gs = ['Garnet', 'Grossular']
mgmaj = ['Garnet', 'MgMajorite']
namaj = ['Garnet', 'NaMajorite']

qz = ['Quartz', 'Quartz']

ky = ['Kyanite', 'Kyanite']

###

fo = ['Olivine', 'Forsterite']
fa = ['Olivine', 'Fayalite']
mgsp = ['Spinel', 'MgSpinel']
hc = ['Spinel', 'Hercynite']
neph = ['Nepheline', 'Nepheline']
mag = ['Ferropericlase','Magnetite']


# plag + opx => gt + cpx
An_CatsQz = ([an], [cats, qz]) # Paria et al (paragraph 2, cites Herzberg 1978) GR67b
AnEn_DiMgtsQz = ([an, en], [di, mgts, qz]) # Green & Ringwood 1967 (a) - GR67a
AnEn_PypDiQz = ([an, en], [pyp, di, qz]) # GR67c, Paria et al (B)
AnEn_PypGsQz = ([an, en], [pyp, gs, qz]) # GR67d
AnFs_AlmHedQz = ([an, fs], [alm, hed, qz]) # Paria et al (A) GR67c
AnFs_AlmGsQz = ([an, fs], [alm, gs, qz]) # GR67d

EnMgts_Pyp = ([en, mgts], [pyp]) # Paria et al (paragraph 2, cites RWood 1974) GR67e, Wood & Banno 1973 (1) Brey et al 1986 (A)
EnCats_GsPyp = ([en, cats], [gs, pyp]) # GR67f
Ab_JdQz = ([ab], [jd, qz]) # Paria et al (paragraph 2, Reinsch 1977 cited in Newton 1983) GR67g

# Gt = Cpx + Qz + Ky
GsPypQz_DiKy = ([gs, pyp, qz], [di, ky]) # Green 1967 (xi) Thermobarometric Methodologies p.246
GsAlmQz_HedKy = ([gs, alm, qz], [hed, ky]) # Green 1967 (xi)
# formation of kyanite
# Plg + Gt = Cpx + ky
An_KyGsQz = ([an], [ky, gs, qz]) # Koziol et al 1988 (paragraph 1) GR67o, Eckert et al. 1991 A

# For AGU2 - allows breakdown of plag without Opx present
# Opx + Tsch = Grt
FsMgts_PypAlm = ([fs, mgts], [pyp, alm]) # analagous to en + mgts => pyp
FsCats_GsPypAlm = ([fs, cats], [gs, pyp, alm]) # analogous to en + cats => gs + pyp
# Cpx + Plg = Grt + Qtz
AnMgts_PypGsQz = ([an, mgts], [pyp, gs, qz]) # MgTs is technically Opx
DiAn_GsPypQz = ([di, an], [gs, pyp, qz]) # Eckert et al. 1991 B

HedAn_GsAlmQz = ([hed, an], [gs, alm, qz])
CatsAn_GsQz = ([cats, an], [gs, qz])

# AGU 3 - bad?
# too much CPX produced at high pressures?
HedCats_GsAlm = ([hed, cats], [gs, alm])
DiMgts_GsPyp = ([di, mgts], [gs, pyp])
HedMgts_GsAlmPyp = ([hed, mgts], [gs, alm, pyp])

# AGU 4
CatsQz_GsKy = ([cats, qz], [gs, ky])
MgtsQz_PypKy = ([mgts, qz], [pyp, ky])

# AGU 5 = undo AGU 3

# AGU 6-8 - didn't seem to work

# AGU 9
EnKy_PypQz = ([en, ky], [pyp, qz]) # Green 67 (x)
FsKy_AlmQz = ([fs, ky], [alm, qz]) # Green 67 (x)

# AGU10 - Mg/Fe exchange - very good
PypHed_AlmDi = ([pyp, hed],[alm, di]) # Ellis & Green 1979 (1)
# "weird" phases
En_Cen = ([en],[cen])
Odi_Di = ([odi],[di])
Mgmaj_En = ([mgmaj],[en])
Namaj_Jd = ([namaj],[jd])

# Added for AGU 11
CatsPyp_GsMgts = ([cats, pyp], [gs, mgts]) # cation exchange with Tschermaks
DiPyp_EnGs = ([di, pyp], [en, gs]) # Brey et al., 1986 (C)
DiCats_GsPyp = ([di, cats], [gs, pyp]) # Brey at al. 1986 (B)

# UNUSED

# Olivine
FoQz_En = ([fo, qz], [en])
FaQz_Fs = ([fa, qz], [fs])
FoEn_Cen = ([fo, en], [cen])
FoEnAn_HedCats = ([fo, en, an], [hed, cats])
FoAn_DiCatsMgts = ([fo, an], [di, cats, mgts])
FoAn_DiCatsEn = ([fo, an], [di, cats, en])
FoAn_DiEnMgts = ([fo, an], [di, en, mgts])
FoAn_CatsEnMgts = ([fo, an], [cats, en, mgts])
FoAn_DiEnMgsp = ([fo, an], [di, en, mgsp])
FoGs_MgspEnDi = ([gs, fo], [mgsp, en, di])
FoPyp_MgspEnDi = ([pyp, fo], [mgsp, en, di])
FoAn_GsPyp = ([fo, an], [gs, pyp])
FaAn_GsAlm = ([fa, an], [gs, alm])

# Spinel
MgspAnEn_GsPyp = ([mgsp, an, en], [gs, pyp])
HcAnFs_GsAlm = ([hc, an, fs], [gs, alm])
AbEnMgsp_JdPyp = ([ab, en, mgsp], [jd, pyp])
AbFsHc_JdAlm = ([ab, fs, hc], [jd, alm])

# ClinoEnstatites
# These are analagous to those above involving MgSiO3 in enstatite
CenMgts_Pyp = ([cen, mgts], [pyp]) # Paria et al (paragraph 2, cites RWood 1974) GR67e
CenCats_GsPyp = ([cen, cats], [gs, pyp]) # GR67f
AnCen_DiMgtsQz = ([an, cen], [di, mgts, qz]) # Green & Ringwood 1967 (a) - GR67a
AnCen_PypDiQz = ([an, cen], [pyp, di, qz]) # GR67c, Paria et al (B)
AnCen_PypGsQz = ([an, cen], [pyp, gs, qz]) # GR67d


agu1 = [ 
    # Simple, based on Green and Ringwood 1967
    An_CatsQz, 
    AnEn_DiMgtsQz, 
    AnEn_PypDiQz, 
    AnEn_PypGsQz,
    AnFs_AlmHedQz,
    AnFs_AlmGsQz,
    EnMgts_Pyp,
    EnCats_GsPyp,
    Ab_JdQz,
]

agu2 = [ 
    # Also based on GR67
    # adds ferrosilite endmember to Opx+Tshermaks => Grt
    FsMgts_PypAlm, 
    FsCats_GsPypAlm, 
    
    # adds Mgts endmember to An+Opx => Grt
    AnMgts_PypGsQz,
    
    # adds a new class of reaction
    # An+Cpx => Grt
    DiAn_GsPypQz, 
    HedAn_GsAlmQz,
    CatsAn_GsQz,
]
    
agu3 = [
    # adds: Cpx + Tshermaks = Garnet
    DiCats_GsPyp,
    HedCats_GsAlm,
    DiMgts_GsPyp,
    HedMgts_GsAlmPyp,

    # These ended up looking wrong, too much CPX produced at high pressures
    # could be other reactions were missing
]

agu4a =  [
    # formation of kyanite
    # An = Ky + Gs + Qz - classic reaction
    An_KyGsQz,
    # Cpx + Ky = Grt + Qz
    GsPypQz_DiKy,
    GsAlmQz_HedKy,
]
agu4b =[
    # related to agu4a
    # Tschermaks + Qz = Grt + Kyanite
    CatsQz_GsKy,
    MgtsQz_PypKy
]

agu6 = [
    # clino/ortho exchange
    En_Cen,
    Odi_Di,

    # important  
]

agu7 = [
    # AGU 6 caused garnet to go into clinoenstatite etc...
    CenMgts_Pyp,
    CenCats_GsPyp,
    AnCen_DiMgtsQz,
    AnCen_PypDiQz,
    AnCen_PypGsQz,

    # these are probably redundant
]

agu9 = [
    # Based on Green '67
    # an+Opx = Cpx + Ky
    # ([an,en],[di,ky]), # redundant, Green 67
    # ([an,fs],[hed,ky]), # redundant, Green 67
    EnKy_PypQz,
    FsKy_AlmQz, 
]

agu10 = [
    # Fe-Mg exchange Cpx-Grt
    # PypHed_AlmDi,
    # majorites
    Mgmaj_En,
    Namaj_Jd,
]

agu11 = [
    CatsPyp_GsMgts, # cation exchange with Tschermaks
    DiPyp_EnGs,
    DiCats_GsPyp, #NOTE: this was in AGU3, which produced too much cpx at the time
]