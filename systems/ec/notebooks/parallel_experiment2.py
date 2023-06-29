from mcm.tcg import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool
import multiprocessing as mp
import importlib
from from_perplex import *
from scipy.integrate import solve_ivp

yr = 3.154e7
kyr = 1e3*yr
Myr = 1e6*yr
s = 1
mm = 1e-3
km = 1e3
g = 1e-3
cm = 1e-2

### ------------ INPUTS -------------------
reference= 'parallel_experiment'
composition = 'hacker_2015_md_xenolith'
rxn_name = 'eclogitization_agu17_stx21_rx'

# only phases greater than this fraction will be plotted
phasetol = 1.e-5 # 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9
max_steps = 3e4

# Set initial and final crustal thickness (m)
moho_i = 30.*km # typical crust
moho_f = 70.*km # Altiplano, Tibet

# Calc descent rate of Moho
shortening_rate = 3.0 *mm/yr # m/s

L0 = 100.*km # lithospheric thickness - m
z0 = moho_i # 
descent_rate = z0/L0 * shortening_rate # m/s

# choose a time step
dt = 100.*kyr # sec
max_t = 100.*Myr
times = np.arange(0,max_t+dt,dt)

depth_m = z0 + times*descent_rate
depth_m[depth_m > 70.*km] = 70.* km

# Calc pressure, temperature from depths
# TODO: use a steady-state geotherm
P_range = 0.027 * depth_m/1e3 # assumes 0.027 Gpa/km
T_moho_i = 800. + 273.15
T_moho_f = 850. + 273.15
T_range = T_moho_i + (depth_m - depth_m[0])/(depth_m[-1] - depth_m[0])*(T_moho_f - T_moho_i)
print(P_range)
show_equilibrium = False
processes =  mp.cpu_count()-1


# initial temperature and pressure in K and bars
T0 = T_moho_i
P0 = P_range[0] # Gpa

# ------------------------------------------

# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-r", "--rxn_name")

    args = parser.parse_args()

    if args.composition is not None:
        composition = args.composition
    if args.rxn_name is not None:
        rxn_name = args.rxn_name

#====================================================

rxn = get_reaction(rxn_name)
phase_names, endmember_names = get_names(rxn)

rxn.set_parameter('T0', T_range[0]) # K; default 2000 K
mi0, Xik0, phii0, Cik0 = get_point_composition(rxn, composition)

ipyrolite = get_rho_interpolator('xu_2008_pyrolite')
iharzburgite = get_rho_interpolator('xu_2008_harzburgite')
rho_pyrolite=ipyrolite((T_range, P_range))
rho_harzburgite=iharzburgite((T_range,P_range))

_Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
_mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0

I = len(rxn.phases())
Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(Kis)

# Equilibrate model at initial (T0, P0)
ode = ScipyPDReactiveODE(rxn)
ode.solve(T_range[0],GPa2Bar(P_range[0]),_mi0,_Cik0,1,Da=1e6,eps=eps)
rho0 = ode.final_rho()*100 # kg/m3

scale= {'T':T0, 'P':P0*1e4, 'rho':rho0, 'h':(moho_f-moho_i)}

# Damkoehler number
Das = [1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]

# NOTE
# reaction rates are multiplied by Damkholer number
# this implies that they are scaled so that mass flux rate = 1
# this defines a time scale?

# Reaction rates based on Hetenyi et al. 2017
# these are based on available surface area of 3/r where r is grain radius
#reaction_rate_gcm =  [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3] # g/cm2/yr

                            #   kg/m3  m/s       
reaction_rate_per_surface = [da*rho0*descent_rate for da in Das]# kg/m2/s
reaction_rate_per_surface_gcm = [r/10*yr for r in reaction_rate_per_surface] # g/cm2/yr

_Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase

def rhs(t,u,rxn,scale,Da,verbose=False):

    # Extract variables
    mi = u[:I]
    Cik = u[I:I+K]
    T = u[I+K] # 1 to 1.something
    P = u[I+K+1] # 1 to 2.something

    # scale Temperature and pressure
    Ts = scale['T']*T
    Ps = scale['P']*P
    print(Ts,Ps)
    C = rxn.zero_C()
    # reshape C
    Kis = 0
    for i,Ki in enumerate(_Kis):
        C[i] = Cik[Kis:Kis+Ki]
        Kis = Kis+Ki
    # regularize C
    Cs = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    Cs = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    # phase properties
    rho = np.array(rxn.rho(Ts, Ps, Cs))
    #alpha = np.array(rxn.alpha(Ts, Ps, C))
    #cp = np.array(rxn.Cp(Ts, Ps, C))
    #s = np.array(rxn.s(Ts, Ps, C))
    v = np.sum(mi/rho)

    # mean properties
    #rho_bar = 1./v
    #alpha_bar = v.dot(alpha)*rho_bar
    #cp_bar = mi.dot(cp)
    #s_bar = mi.dot(s)

    #A = scale['P']/scale['rho']
    
    # regularize m
    mis = np.asarray(mi)
    mis = mi + eps

    Gammai = np.asarray(rxn.Gamma_i(Ts,Ps,Cs,mi))
    gamma_ik = rxn.Gamma_ik(Ts,Ps,Cs,mi)
    Gammaik = np.zeros(K)
    sKi = 0
    for i in range(I):
        for k in range(_Kis[i]):
            Gammaik[sKi+k] = gamma_ik[i][k]
        sKi += _Kis[i]
    
    du = np.zeros(I+K)
    sKi = 0
    for i in range(I):
        du[i] = Da*rho0*Gammai[i]*v
        for k in range(_Kis[i]):
            GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
            du[I+sKi+k] = Da*rho0*GikcGi*v/mis[i]
        sKi += _Kis[i]
    
    dTdz = (T_moho_f - T_moho_i)/scale['T']
    dPdz = 1.8e4/scale['P'] # 0.027*1e4*1e3*(moho_f-moho_i)/scale['h'] # 0.027 GPa/km
    print(dTdz,dPdz)

    #if verbose:
        #print(_T/cp_bar, A*alpha_bar, s.dot(Gammai))
    du = np.empty(u.shape)

    du[I+K:] = np.array([dTdz, dPdz])
    return du

# function to run in parallel
def run_experiment(Da):
    rxn = get_reaction(rxn_name)
    rxn.set_parameter('T0',T0)
    rho_final = np.zeros(T_range.shape)
    phases_final = ['' for i,_ in enumerate(depth_m)]
    phi_final = np.empty(T_range.shape+(I,))
    mi_final = np.zeros(T_range.shape+(I,))
    Cik_final = np.empty(T_range.shape+(K,))
    Xik_final = np.empty(T_range.shape+(K,))

    # Equilibrate model at initial (T0, P0)
    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T_range[0],GPa2Bar(P_range[0]),_mi0,_Cik0,1,Da=1e5,eps=eps)
    
    u0=np.empty(I+K+2)
    u0[:I] = _mi0
    u0[I:I+K] = _Cik0
    u0[I+K:] = np.array([1., 1.]) # T/T0, P/P0
    
    Pmax = (0.027 * moho_f/1e3)/P0 # assumes 0.027 Gpa per km
    event = lambda t, y, rxn, scale,Da: Pmax - P0*y[-1] # when pressure = Pmax
    event.terminal = True
    sol = solve_ivp(rhs,[0,1],u0,args=(rxn,scale,Da),dense_output=True,method='Radau',rtol=rtol,atol=atol, events=event)
    
    T = sol.y[-2][-1]*scale['T']
    P = sol.y[-1][-1]*scale['P']
    print('{} P_end = {} bar.  T_end = {} K. Used {} steps: '.format(sol.message, P,T,len(sol.t)))

    mi  = sol.y[:I,-1]
    Cik = sol.y[I:I+K,-1]

    C  = ode.reshapeC(Cik)
    Cs = ode.regularizeC(C)

    rhoi = np.asarray(rxn.rho(T, P, Cs))

    v    = sum(mi/rhoi)
    rho   = 1/v

    Cik = sol.y[I:I+K,-1]
    mi = sol.y[:I,-1] # -1 = final time step
    mi_reg = ode.regularizem(mi)

    C = ode.reshapeC(Cik)
    Xik = rxn.C_to_X(C)
    Cs = ode.regularizeC(C)
    Cik_reg = [c for arr in Cs for c in arr]
    rhoi = rxn.rho(T, GPa2Bar(P), Cs)
    v = mi/rhoi
    vi  = 1./v.sum()
    phi = vi*mi/rhoi
    rho_final[i]=rho
    phases_final[i]=phases 
    phi_final[i]=phi 
    mi_final[i]=mi_reg
    Cik_final[i]=Cik_reg
    Xik_final[i]= [x for arr in Xik for x in arr]

    return rho_final, phases_final, phi_final, mi_final, Cik_final, Xik_final, Da

# run for varying damkhoeler numbers
with Pool(1) as pool:
    # blocks until all finished
    sols = pool.map(run_experiment, [Das[-1]])

outputPath = Path("figs",reference,composition,rxn_name)
outputPath.mkdir(parents=True, exist_ok=True)

depth_km = depth_m/1e3

# plot 
num_subplots = 2 + len(phase_names) + 2
fig = plt.figure(figsize=(3*num_subplots,12))
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
ax = fig.add_subplot(1,num_subplots,(1,2))
ax.invert_yaxis()

n = len(Das)
colormap = plt.cm.get_cmap('Greys')
ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))
for rho, *rest in sols:
    plt.plot(rho/10, depth_km)

plt.plot(rho_pyrolite/10,depth_km, 'r:')
plt.plot(rho_harzburgite/10,depth_km, 'm:')

plt.legend(['$r_0 = ${:.1e} g/cm$^2$/yr'.format(r) for r in reaction_rate_per_surface_gcm] + (['equil.', 'pyrolite', 'harzburgite'] if show_equilibrium else  ['pyrolite', 'harzburgite']), loc="upper right")

plt.ylabel('Depth (km)')
plt.xlabel('Density')
plt.xlim([2.8, 3.8])

# plot all phases

cmaps = [plt.cm.get_cmap(name) for name in ['Blues', 'YlOrBr', 'Greens','Reds','Purples','copper_r']]
for i,phase in enumerate(phase_names):
    ax = fig.add_subplot(1,num_subplots,3+i)
    ax.invert_yaxis()
    cmap = cmaps[i]
    ax.set_prop_cycle(plt.cycler('color', cmap(np.linspace(0.2, 1, n))))
    ax.set_xlim([0., 80.])
    ax.set_ylabel(None)
    ax.set_xlabel("{} (wt%)".format(phase))
    ax.set_xticks(np.arange(0,110,10))
    ax.set_xticklabels([None, None, 20, None, 40, None, 60, None, 80, None, None])
    ax.set_yticklabels([])
    for rho, phases, phi, mi, Cik, Xik, *rest in sols:
        # plot garnet mass fraction
        plt.plot(mi[:,i]*100, depth_km)
    
# Plot plagioclase composition
ax = fig.add_subplot(1,num_subplots,num_subplots-1)
ax.invert_yaxis()
ax.set_xlim([0., 1.0])
ax.set_ylabel(None)
ax.set_yticklabels([])
ax.set_xlabel("$X_{\mathrm{An}}$")
colormap = plt.cm.get_cmap('Greys')
ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))

for rho, phases, phi, mi, Cik, Xik, *rest in sols:
    Xan = np.asarray(Xik[:,10])
    plt.plot(Xan, depth_km)

# Plot Cpx composition
ax = fig.add_subplot(1,num_subplots,num_subplots)
ax.invert_yaxis()
ax.set_xlim([0., 1.0])
ax.set_ylabel(None)
ax.set_yticklabels([])
ax.set_xlabel("$X_{\mathrm{Jd}}$")
colormap = plt.cm.get_cmap('Greys')
ax.set_prop_cycle(plt.cycler('color', colormap(np.linspace(0.2, 1, n))))

for rho, phases, phi, mi, Cik, Xik, *rest in sols:
    Xjd = np.asarray(Xik[:,4])
    plt.plot(Xjd, depth_km)

fig.suptitle('$R_m=${:.1f} km/Myr, $T_m=${:.0f}-{:.0f} °C, $d_m=${:.0f}-{:.0f} km, {}'.format(descent_rate*1e3,T_moho_i-273.15,T_moho_f-273.15,moho_i/1e3,moho_f/1e3,composition),y=0.9)
plt.savefig(Path(outputPath,"{}.{}".format("results", "pdf")))
plt.savefig(Path(outputPath,"{}.{}".format("results", "png")))
