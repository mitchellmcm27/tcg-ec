# %%

import numpy as np
import matplotlib.pyplot as plt
from fd1d_heat_implicit import *

def rhs():
    return

# From Turcotte & Schubert
# heat production decreases exp from surface
# Math: H = H_0 e^{-y/h_r}


# Subbing into the heat eqn
# Math: 0 = k\frac{d^2T}{dy^2} + \rho H_0 e^{-y/h_r}

# Assume a mantle heat flux q_m

# Integrate equation:
# Math: c_1 = k\frac{dT}{dy} - \rho H_0 h_r e^{-y/h_r} = -q - \rho H_0h_re^{-y/h_r}

# By boundary conditions, c1 = q_m

# Heat flux at any depth
# Math: q = -\left( q_m + \rho H_0h_re^{-y/h_r} \right)

# Surface heat flow:
# Math: q_0 = q_m + \rho h_r H_0

# Sierra Nevada: qm = 17 mW/m2, hr = 10 km
# Eastern US: qm = 33, hr = 7.5
# Norway & Sweden: qm = 22, hr = 7.2
# E Canada Shield: qm = 30.5, hr = 7.1

# hr is usually near 10 km
# qm consistent with mean basal heating of cont'l lithosphere of 28 mW/m2

# Assume T at surface = T0
# Math: T = T_0 + \frac{q_my}{k} + \frac{\rho H_0 h_r^2}{k}(1-e^{-y/h_r})
# or
# Math: T = T_0 + \frac{q_my}{k} + \frac{(q_0 - q_m)h_r}{k}(1-e^{-y/h_r})

def geotherm(depths, thickening):


    T0 = 283 # K, = 10 C
    q0 = 56.6e-3 # W/m2
    qm = 40.e-3/thickening # W/m2
    k = 3.35 # W/m/K
    rhoH0 = 1.5e-6 # W/m3

    # mantle heat prod
    rhoHc = 0.1e-6
    hc = 100.e3*thickening

    hr = 10.e3
    #y = np.asarray(depths)
    #T = T0 + qm*y/k + (rhoH0 * hr)*hr/k*(1-np.exp(-y/hr))
    T = np.zeros(depths.shape)
    for i,y in enumerate(depths):
        if(y<=hc):
            T[i] = T0 + rhoH0*hr*hr/k*(1-np.exp(-y/hr)) - rhoHc*y*y/2/k + (qm + rhoHc*hc)/k*y
        else:
            T[i] = T0 + rhoHc*hc*hc/2/k + rhoH0*hr*hr/k + qm*y/k
    return T

if __name__ == "__main__":

    L0 = 60.e3 # lithosphere thickness, m
    z0 = 30.e3 # moho depth
    thickening = 2.333333 # amount to thicken
    hr0 = 10.e3 # e folding decay depth for radiogenic heat prod, km

    k = 3.35 # thermal conductivity
   
    T_surf = 283. # surface T
    T_lab = 1603. # base lithosphere T
    rhoH0 = 1.9e-6 # W/m3, heat prod at surface
    rhoHc = 0.e-6 # W/m3, heat prod in mantle

    fig = plt.figure()
    ax1 = plt.gca()
    ax1.set_ylabel("Depth (km)")
    ax1.set_xlabel("T (°C)")
    ax2 = plt.gca().twiny()
    ax2.set_xlabel("$\\rho H$ (W/m$^3$)", color="red")
    ax2.tick_params(colors="red",which="both")
    ax2.spines['top'].set_color('red') 

    nsteps = 10
    for (time,L) in enumerate(np.linspace(L0, L0*thickening, nsteps)):
        moho = z0*L/L0
        hr = hr0*L/L0 # initially 10 km
 
        #
        #  X_NUM is the number of equally spaced nodes to use between 0 and 1.
        #
        x_num = 201
        x_min = 0.
        x_max = L
        dx = (x_max - x_min)/(x_num - 1)
        x = np.linspace(x_min, x_max, x_num)
        #
        #  T_NUM is the number of equally spaced time points between 0 and 10.0.
        #
        t_num = 101
        t_min = 0.0
        t_max = 1e6 * 3.156e7
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
            H = np.zeros(X.shape)
            for i,x in enumerate(X):
                if(x > moho):
                    H[i] = rhoH0*np.exp(-x/hr) # uniform mantle H
                else:
                    H[i] = rhoH0*np.exp(-x/hr) # decaying H within crust
            return H
        
        def bc_fun(x_num,X,t,u):
            u[0] = T_surf
            u[x_num-1] = T_lab
            return u
        
        for j in range(0, t_num):
            if (j == 0 and time==0):
                u = np.linspace(T_surf,T_lab,x_num) # initial condition
            else:
                u = fd1d_heat_implicit(a,x_num,x,t[j-1],dt,cfl,rhs_fun,bc_fun,u)
            Tmoho = np.interp(moho,x,u)
        #
        #  Plot X and T versus H.
        #
        ax1.plot(u[x<moho]-273.,x[x<moho]/1.e3, 'k', linewidth=1,alpha=time/nsteps/2 + 0.4)
        ax1.plot(u[x>=moho]-273.,x[x>=moho]/1.e3, 'k', linewidth=1,alpha=time/nsteps/2 + 0.1)
        ax1.plot(Tmoho-273.,moho/1e3,'ko')
        ax2.plot(rhoH0*np.exp(-x/hr),x/1.e3, 'r', linewidth=1, alpha=time/nsteps/2 + 0.1)
    ax1.invert_yaxis()
    plt.title("$L_0=${:n} km, $z_0=${:n} km, $h_{{r0}}=${:n} km\n$T_s=${:n} °C, $T_L=${:n} °C, $\\rho H_s=${:g} μW/m$^3$, $k=${}".format(L0/1e3,z0/1e3,hr0/1e3,T_surf-273, T_lab-273, rhoH0*1e6,k),
              fontsize=8)
    plt.show()

# %%
