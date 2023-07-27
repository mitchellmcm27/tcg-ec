import numpy as np
import matplotlib.pyplot as plt
from fd1d_heat_implicit import *

def geotherm_steady(m, # between 0 and 1
                    L, # length scale
                    thickening = 1.0, # gt 1
                    Ts=283.15,
                    Tlab=1603.15,
                    k=3.35,
                    A=1.9e-6,
                    hr0=10.e3):
    
    hr = hr0 * thickening
    b1 = Ts
    b2 = Tlab

    b3 = hr**2*A/k
    b4 = hr*A*L/k

    a1 = b1 + b3
    a0 = b2 - a1 + b3*np.exp(-L/hr)
    T = a0*m + a1 - b3*np.exp(-L*m/hr)

    dTdx = lambda x: a0 + b4*np.exp(-L*x/hr)
    q0 = k * dTdx(0)/L

    return T, q0

if __name__ == "__main__":
    T, q0 = geotherm_steady(30./60., 60.e3, A=2.5e-6, k=2.25, hr0=14.e3)
    print("Hot: T_moho={:.0f} C, q0={:.0f} mW/m2".format(T-273.15,q0*1e3))
    
    T, q0 = geotherm_steady(30./70., 70.e3, A=2.5e-6, k=2.25, hr0=10.e3)
    print("Transitional: T_moho={:.0f} C, q0={:.0f} mW/m2".format(T-273.15,q0*1e3))
    
    T, q0 = geotherm_steady(30./80., 80.e3,A=2.5e-6,k=2.25,hr0=10.e3)
    print("Cold: T_moho={:.0f} C, q0={:.0f} mW/m2".format(T-273.15,q0*1e3))
