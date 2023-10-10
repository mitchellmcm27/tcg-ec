import numpy as np
import matplotlib.pyplot as plt
from fd1d_heat_implicit import *
from scipy.integrate import solve_bvp

def geotherm_steady(m, # between 0 and 1
                    L, # length scale
                    thickening = 1.0, # gt 1
                    Ts=283.15, # K
                    Tlab=1603.15, #K
                    k=3.35, # W/m/K
                    A=1.9e-6, # W/m3
                    hr0=10.e3, # km
                    method="analytic"):
    hr = hr0 * thickening

    if method=="analytic":
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
    else:
        z = np.linspace(0,1,1000)*L*thickening # [0, 0.001, 0.002, ..., 1]
 
        def f(x,y):
            return np.vstack((y[1], np.full_like(x, -A*np.exp(-x/hr)/k)))
        def bc(ya,yb):
            return np.array([ya[0]-Ts, yb[0]-Tlab, ya[1], yb[1]-0.4])

        # initial guess
        y0 = np.asarray([np.linspace(Ts,Tlab,1000),np.ones(1000)])
        sol = solve_bvp(f, bc, z, y0)
        T = sol.y[0]
        dTdx = sol.y[1]
        q0 = k * dTdx[0]/L
        print(T)
        plt.plot(T,-z/1e3)
        plt.show()
    return T[0], q0

if __name__ == "__main__":
    T, q0 = geotherm_steady(30./60., 60.e3, A=2.5e-6, k=2.25, hr0=14.e3)
    print("Hot: T_moho={:.0f} C, q0={:.0f} mW/m2".format(T-273.15,q0*1e3))
