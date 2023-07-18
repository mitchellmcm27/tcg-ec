import numpy as np
from fd1d_heat_implicit import *

T_surf = 283.15 # surface T, K
T_lab = 1603.15 # base lithosphere T, K
rhoH0 = 1.9e-6 # W/m3, heat prod at surface
k = 3.35
hr0 = 10.e3
def get_geotherm(L0,z0,thickening,t_max):
    L = L0*thickening
    moho = z0*L/L0
    hr = hr0*L/L0

    #
    #  X_NUM is the number of equally spaced nodes to use between 0 and 1.
    #
    x_num = 201
    x_min = 0.
    x_max = L
    x = np.linspace(x_min, x_max, x_num)

    if (moho not in x):
        x = np.insert(x, moho).sort()
        x_num += 1
    moho_idx = np.where(x == moho)[0][0]
    dx = (x_max - x_min)/(x_num - 1)

    #
    #  T_NUM is the number of equally spaced time points between 0 and 10.0.
    #
    t_num = 101
    t_min = 0.0
   
    dt = (t_max - t_min)/(t_num - 1)
    t = np.linspace(t_min, t_max, t_num)
    x_step = (x_max - x_min)/float(x_num - 1)
    t_step = (t_max - t_min)/float(t_num - 1)
    cfl = k*t_step/x_step/x_step

    # system matrix
    a = np.zeros (( x_num, x_num ))
    a[0,0] = 1.0
    for i in range ( 1, x_num - 1 ):
        a[i,i-1] = -cfl
        a[i,i  ] = 1.0 + 2.0 * cfl
        a[i,i+1] = -cfl
    a[x_num-1,x_num-1] = 1.0

    def rhs_fun(x_num, X, t):         
        return rhoH0*np.exp(-X/hr) # decaying H within crust
    
    def bc_fun(x_num,X,t,u):
        u[0] = T_surf
        u[x_num-1] = T_lab
        return u
    
    for j in range(0, t_num):
        if (j == 0):
            u = np.linspace(T_surf,T_lab,x_num) # initial condition
        else:
            u = fd1d_heat_implicit(a,x_num,x,t[j-1],dt,cfl,rhs_fun,bc_fun,u)
    return u[moho_idx]
    
