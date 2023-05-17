# Returns an interpolant for density of various mantle compositions in kg/m3/100
# Rho from equilibrium data computed in Perple_X 
# Compositions from Xu et al., 2008 EPSL

from scipy.interpolate import LinearNDInterpolator
import pandas as pd
import numpy as np

def get_rho_interpolator(filepath):
    df = pd.read_csv(filepath, delimiter='\s+', header=0, names=["T","P","rho"])
    P = df["P"].to_numpy()
    T = df["T"].to_numpy()
    rho = df["rho"].to_numpy()/100
    interp = LinearNDInterpolator((T,P), rho)
    return interp