from mcm.tcg import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import *
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from multiprocessing import Pool, shared_memory
import multiprocessing as mp
from from_perplex import *
from scipy.integrate import solve_ivp
from geotherm_steady import geotherm_steady

yr = 3.154e7
kyr = 1e3*yr
Myr = 1e6*yr
s = 1
mm = 1e-3
km = 1e3
g = 1e-3
cm = 1e-2

### ------------ INPUTS -------------------
reference= "parallel_experiment2"
rxn_name = "eclogitization_agu17_stx21_rx"

# only phases greater than this fraction will be plotted
phasetol = 1.e-5 # 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9
max_steps = 3e4

# Calc descent rate of Moho
shortening_rate = 3.0 *mm/yr # m/s

crustal_rho = 2780.
gravity = 9.81 

tectonic_settings = [
    {
        "setting": "orogen_hot",
        "L0": 60.e3,
        "z0": 30.e3,
        "z1": 70.e3,
        "As": 2.5e-6,
        "hr0": 10.e3,
        "k": 2.25,
    },
    {
        "setting": "orogen_cold",
        "L0": 100.e3,
        "z0": 35.e3,
        "z1": 70.e3,
        "As": 2.0e-6,
        "hr0": 10.e3,
        "k": 2.25,
    },
    {
        "setting":"craton_hot",
        "L0": 125.e3,
        "z0": 35.e3,
        "z1": 70.e3,
        "As": 3.0e-6,
        "hr0": 15.e3,
        "k": 2.25,
    },
    {
        "setting":"craton_cold",
        "L0": 135.e3,
        "z0": 35.e3,
        "z1": 70.e3,
        "As": 2.0e-6,
        "hr0": 10.e3,
        "k": 2.25,
    }
]

# Damkoehler numbers
Das = [1e-3]#, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]#, 1e5]
compositions = ["hacker_2015_md_xenolith", "sammon_2021_deep_crust", "zhang_2022_cd07-2"]

scenarios = []

for setting in tectonic_settings:
    for da in Das:
        for comp in compositions:
            scenario = setting.copy()
            scenario["Da"] = da
            scenario["composition"] = comp
            scenario["Ts"] = 10. + 273.15
            scenario["Tlab"] = 1330. + 273.15
            scenarios.append(scenario)

processes =  mp.cpu_count()-1
# ------------------------------------------

# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-r", "--rxn_name")

    args = parser.parse_args()

    if args.composition is not None:
        compositions = [args.composition]
    if args.rxn_name is not None:
        rxn_name = args.rxn_name

#====================================================

rxn = get_reaction(rxn_name)

phase_names, endmember_names = get_names(rxn)
I = len(rxn.phases())
Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(Kis)

ic_cache = {}
def get_initial_composition(T0,P0,rxn,composition):
    slug = "{:.4f}-{:.4f}-{}-{}".format(T0,P0,rxn.name(),composition)
    if slug in ic_cache:
        print("ic_cache hit")
        return ic_cache[slug]

    # Get equilibrium composition at arbitrary P,T
    mi0, Xik0, phii0, Cik0 = get_point_composition(rxn, composition) 
    _Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
    _mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0

    # Equilibrate model at initial (T0, P0)
    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T0, P0,_mi0,_Cik0,1,Da=1e6,eps=eps)
    
    rho0 = ode.final_rho()*100 # kg/m3
    Cik0 = ode.sol.y[ode.I:ode.I+ode.K,-1] # -1 = final time step
    mi0 = ode.sol.y[:ode.I,-1]

    ic_cache[slug] = (Cik0, mi0, rho0)
    return Cik0, mi0, rho0

# initial temperature and pressure in K and bars
def get_initial_state(rxn,scenario):

    composition = scenario['composition']
    L0 = scenario["L0"]
    z0 = scenario["z0"]
    As = scenario["As"]
    hr0 = scenario["hr0"]
    k = scenario["k"]
    Ts = scenario["Ts"]
    Tlab = scenario["Tlab"]

    # Temperature
    T0, qs0 = geotherm_steady(
        z0/L0, # between 0 and 1
        L0,
        thickening=1.0, # gt 1
        Ts=Ts,
        Tlab=Tlab,
        k=k,
        A=As,
        hr0=hr0
    )
    # Pressure
    P0 = crustal_rho * gravity * z0/1e5 # bar

    Cik0, mi0, rho0 = get_initial_composition(T0,P0,rxn,composition)

    return T0, P0, Cik0, mi0, rho0

ipyrolite = get_rho_interpolator("xu_2008_pyrolite")
iharzburgite = get_rho_interpolator("xu_2008_harzburgite")

_Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase

def rhs(t,u,rxn,scale,Da,L0,z0,As,hr0,conductivity,T_surf,Tlab):

    # Extract variables
    mi = u[:I]
    Cik = u[I:I+K]
    T = u[I+K]
    P = u[I+K+1] 
 
    # scale Temperature and pressure
    Ts = scale["T"]*T
    Ps = scale["P"]*P
    dz = scale["h"] # m per unit time
    rho0 = scale["rho"]

    C = rxn.zero_C()
    
    # reshape C
    Kis = 0
    for i,Ki in enumerate(_Kis):
        C[i] = Cik[Kis:Kis+Ki]
        Kis = Kis+Ki

    # regularize C
    Cs = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    Cs = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    rhoi = np.array(rxn.rho(Ts, Ps, Cs))
    v = np.sum(mi/rhoi)
    
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
    
    du = np.zeros(u.shape)
    sKi = 0
    for i in range(I):
        du[i] = Da*rho0*Gammai[i]*v
        for k in range(_Kis[i]):
            GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
            du[I+sKi+k] = Da*rho0*GikcGi*v/mis[i]
        sKi += _Kis[i]
    
  
    # linear temperature
    # dT_linear = (T_moho_f - T_moho_i)/scale["T"]

    dP = crustal_rho * gravity / 1e5 * dz # bar 
    dP = dP/scale["P"]

    shortening = 1 + ((dz+z0)/z0 - 1)*t
    T_steady, q_s = geotherm_steady(z0/L0,
                                    L0*shortening,
                                    shortening,
                                    Ts=T_surf,
                                    Tlab=Tlab,
                                    k=conductivity,
                                    A=As,
                                    hr0=hr0)
    dT = (T_steady-Ts)*dz/scale["T"]

    du[I+K:] = np.array([dT, dP])
    
    return du

def setup_ics(scenario):
    rxn = get_reaction(rxn_name)
    T0, P0, Cik0, mi0, rho0 = get_initial_state(rxn,scenario)
    scenario["T0"] = T0
    scenario["P0"] = P0
    scenario["Cik0"] = Cik0
    scenario["mi0"] = mi0
    scenario["rho0"] = rho0
    return scenario

# function to run in parallel
def run_experiment(scenario):

    # print("Running ", scenario, "...\n")
    Da = scenario["Da"]
    L0 = scenario["L0"]
    z0 = scenario["z0"]
    As = scenario["As"]
    hr0 = scenario["hr0"]
    z1 = scenario["z1"]
    T0 = scenario["T0"]
    P0 = scenario["P0"] 
    Cik0 = scenario["Cik0"]
    mi0 = scenario["mi0"]
    rho0 = scenario["rho0"]
    k = scenario["k"]
    Ts = scenario["Ts"]
    Tlab = scenario["Tlab"]

    rxn = get_reaction(rxn_name)
    rxn.set_parameter("T0",T0)

    scale= {"T":T0, "P":P0, "rho":rho0, "h":(z1-z0)}
    print(scale)

    u0=np.empty(I+K+2)
    u0[:I] = mi0
    u0[I:I+K] = Cik0
    u0[I+K:] = np.array([1., 1.]) # T/T0, P/P0
    
    #Pmax = crustal_rho * gravity * z1/1e5 # bar
    #print(Pmax)
    #event = lambda t,y,rxn,scale,Da: Pmax - y[-1]*scale["P"] # when pressure = Pmax
    #event.terminal = True
    args = (rxn,scale,Da,L0,z0,As,hr0,k,Ts,Tlab)
    sol = solve_ivp(rhs, [0,1], u0, args=args, dense_output=True, method="BDF", rtol=rtol, atol=atol, events=None)
    
    t = np.linspace(0,1,1000)
    y = sol.sol(t)

    T = y[-2]*scale["T"] # K
    P = y[-1]*scale["P"] # bar

    mi_times  = y[:I].T
    Cik_times = y[I:I+K].T

    # reshape C
    def reshape_C(Cik):
        C = rxn.zero_C()
        k = 0
        for i,Ki in enumerate(Kis):
            C[i] = Cik[k:k+Ki]
            k = k+Ki
        return C
    
    Cs_times = [reshape_C(Cik) for Cik in Cik_times]

    rho = [1/sum(mi_times[t]/rxn.rho(T[t], P[t], Cs))/10 for t,Cs in enumerate(Cs_times)]

    print("{} P_end = {:.2f} Gpa. T_end = {:.2f} K. Used {:n} steps: ".format(sol.message,P[-1]/1e4,T[-1],len(sol.t)))

    depth_m = (P*1e5) / gravity / crustal_rho
    
    scenario["T"] = T # K
    scenario["P"] = P # bar
    scenario["rho"] = rho # g/cm3
    scenario["mi"] = mi_times
    scenario["Cik"] = Cik_times
    scenario["Xik"] = np.asarray([rxn.C_to_X(c) for c in Cs_times], dtype="object")
    scenario["z"] = depth_m

    return scenario

print("Preparing to run {} scenarios".format(len(scenarios)))

scenarios = [setup_ics(s) for s in scenarios]

# run for varying damkhoeler numbers
with Pool(processes) as pool:
    # blocks until all finished
    outs = pool.map(run_experiment, scenarios)

for composition in compositions:
    for tectonic_setting in tectonic_settings:

        setting = tectonic_setting["setting"]

        # same setting, same composition
        # different Da
        
        outs_c = sorted([out for out in outs if (out["composition"] == composition and out["setting"] == setting)], key=lambda out: out["Da"])

        base = outs_c[0]
        T = base["T"]
        P = base["P"]
        rho0 = base["rho0"]
        L0 = base["L0"]
        z0 = base["z0"]
        z1 = base["z1"]

        _Das = [o["Da"] for o in outs_c]
        
        descent_rate = z0/L0 * shortening_rate # m/s
        max_t = (z1-z0)/descent_rate
        r0 = [da*rho0/max_t for da in _Das] # Gamma0 (kg/m3/s)
        dg = 3.*mm # grain size
        reaction_rate_per_surface = [r*dg for r in r0] # r0 (kg/m2/s), assumes 1mm grains
        reaction_rate_per_surface_gcm = [r/10*yr for r in reaction_rate_per_surface] # g/cm2/yr

        #print(reaction_rate_per_surface_gcm)
        output_path = Path("figs",reference,rxn_name,setting,composition)
        output_path.mkdir(parents=True, exist_ok=True)

        # plot 
        num_subplots = 3 + len(phase_names) + 2
        subplot_mosaic = [["rho","rho"] + phase_names + ["An", "Jd", "T"]]

        fig = plt.figure(figsize=(3*num_subplots,12))
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        axes = fig.subplot_mosaic(subplot_mosaic)
        [ax.invert_yaxis() for label,ax in axes.items()]

        num_lines = len(_Das)
   
        greys = plt.cm.get_cmap("Greys")
        axes["rho"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))
        axes["An"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))
        axes["Jd"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))

        for i, obj in enumerate(outs_c):
            ax = axes["rho"]
            ax.plot(obj["rho"], obj["z"])

        rho_pyrolite=ipyrolite((T, P/1e4))
        rho_harzburgite=iharzburgite((T,P/1e4))

        ax.plot(rho_pyrolite/10, outs_c[0]["z"], "r:")
        ax.plot(rho_harzburgite/10, outs_c[0]["z"], "m:")
        ax.legend(["$r_0 = ${:.1e} g/cm$^2$/yr".format(r) for r in reaction_rate_per_surface_gcm] +  (["pyrolite", "harzburgite"]), loc="upper right")

        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Density")
        ax.set_xlim([2.8, 3.8])

        # Plot temperature profile
        ax = axes["T"]
        ax.plot(T-273.15, base["z"]/1e3, linewidth=2)
        ax.set_xlabel("T (°C)")
        ax2 = ax.twinx()
        ax2.plot(T-273.15, P/1e4, alpha=0)
        ax2.invert_yaxis()
        ax2.set_ylabel("P (GPa)")

        # plot all phases
        cmaps = [plt.cm.get_cmap(name) for name in ["Blues", "YlOrBr", "Greens","Reds","Purples","copper_r"]]
        for i,phase in enumerate(phase_names):
            ax = axes[phase]
            cmap = cmaps[i]
            ax.set_prop_cycle(plt.cycler("color", cmap(np.linspace(0.2, 1, num_lines))))
            ax.set_xlim([0., 80.])
            ax.set_ylabel(None)
            ax.set_xlabel("{} (wt%)".format(phase))
            ax.set_xticks(np.arange(0,110,10))
            ax.set_xticklabels([None, None, 20, None, 40, None, 60, None, 80, None, None])
            ax.set_yticklabels([])
            for j, obj in enumerate(outs_c):
                ax.plot(obj["mi"][:,i]*100, obj["z"])

        # Plot plagioclase composition
        ax = axes["An"]
        ax.set_xlim([0., 1.0])
        ax.set_ylabel(None)
        ax.set_yticklabels([])
        ax.set_xlabel("$X_{\mathrm{An}}$")

        ax = axes["Jd"]
        ax.set_xlim([0., 1.0])
        ax.set_ylabel(None)
        ax.set_yticklabels([])
        ax.set_xlabel("$X_{\mathrm{Jd}}$")


        for j, obj in enumerate(outs_c):
            XAn = [x[3][0] for x in obj["Xik"]]
            XJd = [x[0][4] for x in obj["Xik"]]

            axes["An"].plot(XAn, obj["z"])
            axes["Jd"].plot(XJd, obj["z"])

        fig.suptitle("{}, {}, $R_m=${:.1f} km/Myr, $d_g=${}mm".format(composition.capitalize().replace("_"," "), setting, descent_rate/1e3*yr*1e6, dg/mm),y=0.9)
        plt.savefig(Path(output_path,"{}.{}".format("results", "pdf")))
        plt.savefig(Path(output_path,"{}.{}".format("results", "png")))
