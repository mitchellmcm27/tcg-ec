# Returns an interpolant for density of various mantle compositions in kg/m3/100
# Rho from equilibrium data computed in Perple_X 
# Compositions from Xu et al., 2008 EPSL

from scipy.interpolate import LinearNDInterpolator,RegularGridInterpolator
import pandas as pd
import numpy as np
from io import StringIO

def get_rho_interpolator(name):
    filepath = "perple_x/output/{}/{}_2.tab".format(name,name)
    try:
        with open(filepath, 'r') as fp:
            '''
            |6.6.6
            sammon_2021_lower_crust_2.tab
            2
            T(K)
            773.14999999999998
            3.1446540880503151
            160
            P(bar)
            5000.0000000000000
            125.78616352201257
            160
            3
            '''         
            perplex_version = fp.readline()
            perplex_fname = fp.readline()
            n_dims = fp.readline() # ?
            xvarname = fp.readline()
            x0 = float(fp.readline().strip())
            xstep = float(fp.readline().strip())
            nx = int(fp.readline().strip())
            yvarname = fp.readline()
            y0 = float(fp.readline().strip())
            ystep = float(fp.readline().strip())
            ny = int(fp.readline().strip())
            n_cols = fp.readline() # ?
            T_range = np.arange(0,nx)*xstep + x0
            P_range = np.arange(0,ny)*ystep + y0
        df = pd.read_csv(filepath, delimiter='\s+', skiprows=12, header=0, names=["T","P","rho"])
    except Exception as e:
        print(e)
        return None
 
    rho = df["rho"].to_numpy()/100 # convert to 100 kg/m3
    rho_g = np.reshape(rho,(nx,ny), order='F')
    interp = RegularGridInterpolator((T_range, P_range/1e4), rho_g, bounds_error=False, fill_value=np.nan)
    return interp

def get_profile_data(name):
    filepath = "perple_x/output/{}/{}_1.tab".format(name,name)
    try:
        df = pd.read_csv(filepath, delimiter='\s+', skiprows=8, header=0)
        df.fillna(0, inplace=True)
        pl0 = df["Pl"]
        pl1 = 0 if not "Pl.1" in df.columns else df["Pl.1"]
        cpx0 = df["Cpx"]
        cpx1 = 0 if not "Cpx.1" in df.columns else df["Cpx.1"]
        cpx2 = 0 if not "Cpx.2" in df.columns else df["Cpx.2"]
        gt0 = df["Gt"]
        gt1 = 0 if not "Gt.1" in df.columns else df["Gt.1"]
        df["Pl2"] = pl0 + pl1
        df["Cpx3"] = cpx0 + cpx1 + cpx2
        df["Gt2"] = gt0 + gt1
        return df
    except Exception as e:
        print("Error parsing perple_x profile:")
        print(e)
        return None


def get_point_composition(rxn, name):

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
        ms = m_df['wt%']
        mi0 = [0 for p in rxn.phases()]

        for n,p in enumerate(rxn.phases()):
            if p.abbrev()=="cpx":
                mi0[n] = ms.get("Cpx",0.)
            if p.abbrev()=="opx":
                mi0[n] = ms.get("Opx",0.)
            if p.abbrev()=="qtz":
                mi0[n] = ms.get("qtz",0.)
            if p.abbrev()=="plg":
                mi0[n] = ms.get("Pl",0.)
            if p.abbrev()=="gt":
                mi0[n] = ms.get("Gt",0.)
            if p.abbrev()=="ky":
                mi0[n] = ms.get("ky",0.)
            if p.abbrev()=="sp":
                mi0[6] = ms.get("Sp",0.)
            if p.abbrev()=="co":
                mi0[6] = ms.get("Aki",0.)

        mi0 = [m/100 for m in mi0]
        #print(mi0)
        x = x.strip().split("\n")

        Xik0 = rxn.zero_C()
        for i,c in enumerate(Xik0):
            Xik0[i][0] = 1.

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
            elif(phase=="Sp"):
                Xik0[6] = [ems[1],ems[0]]
            elif(phase=="Cpx"):
                Xik0[0] = [ems[1], ems[2], ems[3], ems[4], ems[0]]
            elif(phase=="Opx"):
                Xik0[1] = [ems[1], ems[2], ems[3], ems[0]]
            elif(phase=="Gt"):
                if(len(ems)==5):
                    # 5 endmember garnet
                    Xik0[4] = [ems[4], ems[3], ems[2], ems[1], ems[0]]
                    #print(Xik0)
                    # regularize 3-component garnet
                    g3 = (1-(Xik0[4][0]+Xik0[4][1]+Xik0[4][2]))/3.0
                    Xik0[4][0] += g3
                    Xik0[4][1] += g3
                    Xik0[4][2] += g3
                    Xik0[4][3] = 0.0
                    Xik0[4][4] = 0.0

                elif(len(ems)==3):
                    # 3 endmember garnet
                    Xik0[4] = [ems[2], ems[1], ems[0], 0., 0.,]

        phii0 = None
        Cik0 = None
        return mi0, Xik0, phii0, Cik0
    except Exception as e:
        print("There was an error:")
        print(e)
        return None