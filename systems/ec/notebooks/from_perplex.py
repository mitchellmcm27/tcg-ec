# Returns an interpolant for density of various mantle compositions in kg/m3/100
# Rho from equilibrium data computed in Perple_X 
# Compositions from Xu et al., 2008 EPSL

from scipy.interpolate import LinearNDInterpolator
import pandas as pd
import numpy as np
from io import StringIO

def get_rho_interpolator(name):
    filepath = "perple_x/output/{}/{}_2.tab".format(name,name)
    print(filepath)
    try:
        df = pd.read_csv(filepath, delimiter='\s+', skiprows=12, header=0, names=["T","P","rho"])
    except:
        return None

    P = df["P"].to_numpy()
    T = df["T"].to_numpy()
    rho = df["rho"].to_numpy()/100
    interp = LinearNDInterpolator((T,P), rho)
    return interp

def get_profile_data(name):
    filepath = "perple_x/output/{}/{}_1.tab".format(name,name)
    try:
        df = pd.read_csv(filepath, delimiter='\s+', skiprows=8, header=0)
        df.fillna(0, inplace=True)
        df["Pl2"] = df["Pl"] + df["Pl.1"]
        return df
    except Exception as e:
        print(e)
        return None

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


def get_point_composition(name):

    filepath = "perple_x/output/{}/{}_1.txt".format(name,name)
    print(filepath)
    try:
        x = ""
        m = ""
        with open(filepath) as file:
            copy = False
            for line in file:
                if line.strip().startswith("Phase Compositions"):
                    copy = True
                elif line.strip().startswith("Phase speciation"):
                    copy = False
                elif copy:
                    m += line.strip() + "\n"
        with open(filepath) as file:
            copy = False
            for line in file:
                if line.strip().startswith("Phase speciation"):
                    copy = True
                elif line.strip().startswith("Structural formulae for 688 format solution models"):
                    copy = False
                elif copy:
                    x += line.strip() + "\n"
  
        m = m.strip().split("\n")
        m[0] = m[0].replace(" %", "%")
        m = "\n".join(m)
        m_df = pd.read_csv(StringIO(m), delimiter="\s+", header=0)

        mi0 = [0 for p in phase_names]
        ms = m_df["wt%"]
        mi0[0] = ms.get("Cpx",0.)
        mi0[1] = ms.get("Opx",0.)
        mi0[2] = ms.get("qtz",0.)
        mi0[3] = ms.get("Pl",0.)
        mi0[4] = ms.get("Gt",0.)
        mi0[5] = ms.get("ky",0.)

        mi0 = [m/100 for m in mi0]
        #print(mi0)
        x = x.strip().split("\n")

        Xik0 = [
            [1., 0., 0., 0., 0.], # di, hed, cEn, cats, jd
            [1., 0., 0., 0.], # en, fs, mgts, oDi
            [1.], # quartz
            [1.0, 0.], # an, ab
            [1., 0., 0., 0., 0.], # py, alm, gr, *mgmaj, *namaj
            [1.], # kyanite
        ]

        import re
        for line in x:
            l = re.sub(r'\s+',' ', line.strip())
            [phase, rest] = l.split(" ", 1)
            ems = []
            for vals in rest.split(','):
                [k,v] = vals.strip().split(" ")
                ems.append(float(v))
            if(phase=="Pl"):
                Xik0[3] = [ems[1],ems[0]]
            elif(phase=="Cpx"):
                Xik0[0] = [ems[1], ems[2], ems[3], ems[4], ems[0]]
            elif(phase=="Opx"):
                Xik0[1] = [ems[1], ems[2], ems[3], ems[0]]
            elif(phase=="Gt"):
                Xik0[4] = [ems[4], ems[3], ems[2], ems[1], ems[0]]
        #print(Xik0)
        # regularize 3-component garnet
        g3 = (1-(Xik0[4][0]+Xik0[4][1]+Xik0[4][2]))/3.0
        Xik0[4][0] += g3
        Xik0[4][1] += g3
        Xik0[4][2] += g3
        Xik0[4][3] = 0.0
        Xik0[4][4] = 0.0

        phii0 = None
        Cik0 = None
        return mi0, Xik0, phii0, Cik0
    except Exception as e:
        print("There was an error:")
        print(e)
        return None