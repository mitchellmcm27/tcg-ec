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
    Te = lambda x: a0*x + a1 - b3*np.exp(-L*x/hr)
    T = Te(m)

    dTdx = lambda x: a0 + b4*np.exp(-L*x/hr)
    q0 = k * dTdx(0)/L

    return T, q0

if __name__ == "__main__":
    T, q0 = geotherm_steady(0.5,80.e3,A=2.0e-6)
    print("Hot orogen: T_moho={:.0f} C, q0={:.0f} mW/m2".format(T-273.15,q0*1e3))

    T, q0 = geotherm_steady(.3,100.e3,A=2.0e-6)
    print("Cold orogen: T_moho={:.0f} C, q0={:.0f} mW/m2".format(T-273.15,q0*1e3))

    # Gray & Pysklywec?
    T, q0 = geotherm_steady(0.26667,150.e3,A=2.0e-6)
    print("Craton: T_moho={:.0f} C, q0={:.0f} mW/m2".format(T-273.15,q0*1e3))
