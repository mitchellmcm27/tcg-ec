from mcm.tcg import *
import numpy as np
import numpy.ma as ma
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
import csv

####################
# Unit conversions #
####################

yr = 3.154e7
kyr = 1e3*yr
Myr = 1e6*yr
s = 1
mm = 1e-3
km = 1e3
g = 1e-3
cm = 1e-2

##########
# Inputs #
##########

save_output = False
load_output = True

reference= "parallel_experiment2"
rxn_name = "eclogitization_agu17_stx21_rx"

# only phases greater than this fraction will be plotted
phasetol = 1.e-5 # default 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # default 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9

# Calc descent rate of Moho
shortening_rate = 10.0 *mm/yr # m/s

# thickening rate depends on the 
# width of the region being shortened
# and the thickness of the lithosphere
# here we assume that the thickness is 1/3 the width
# e.g. 100 km lithospehre, 300 km region
# this gives a thickening rate of ~3 km/Myr
# and a crustal thickening rate of ~1 km/Myr
# which roughly matches Wang et al. 2015 Arizaro
thickening_rate = shortening_rate / 3. # m/s

# when crust is ~ 1/3 of lithosphere, this gives a descent rate of approx:
moho_descent_rate = 1.0 * mm/yr
h0 = 50. * km # thicken the crust by 50 km
v0 = moho_descent_rate
t0 = h0 / v0

crustal_rho = 2780.
gravity = 9.81

# multiprocessing
num_processes =  mp.cpu_count()

# allows deterministic PDFs
pdf_metadata = {'CreationDate': None}

# Damkoehler numbers
Das = [1e-2, 3e-2, 1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]

# default end time (scaled) is 1
end_t = 1.

# prefix file path for saving plots
prefix = None

# Compositions
compositions = [
    "hacker_2015_bin_4",
    "hacker_2015_bin_3",
    #"bhowany_2018_hol2a",
    #"zhang_2006_mafic_granulite",
    "sammon_2021_lower_crust_norm",
    "sammon_2021_deep_crust_norm",
    "hacker_2015_bin_2",
    "hacker_2015_md_xenolith_norm",
    "hacker_2015_bin_1",
    #"zhang_2022_cd07-2",
    "mackwell_1998_maryland_diabase_norm"
]

compositions = [
    "si50_norm",
    "si51_norm",
    "si52_norm",
    "si53_norm",
    "si54_norm",
    "si55_norm",
    "si56_norm",
    "si57_norm",
    "si58_norm"
]

tectonic_settings = [
    {
        "setting": "hot-1",
        "L0": 55.e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 2.0e-6,
        "hr0": 13.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },    
    {
        "setting": "hot-2",
        "L0": 60.e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.95e-6,
        "hr0": 12.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "hot-3",
        "L0": 65.5e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.9e-6,
        "hr0": 12.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "transitional-1",
        "L0": 71.5e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.85e-6,
        "hr0": 11.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "transitional-2",
        "L0": 78.e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.8e-6,
        "hr0": 11.0e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "transitional-3",
        "L0": 85.e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.75e-6,
        "hr0": 10.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "cold-1",
        "L0": 92.5e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.7e-6,
        "hr0": 10.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting":"cold-2",
        "L0": 102.5e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.65e-6,
        "hr0": 9.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting":"cold-3",
        "L0": 111.e3,
        "z0": 30.e3,
        "z1": 80.e3,
        "As": 1.6e-6,
        "hr0": 9.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    }
]


#######################
# Parse CLI arguments #
#######################

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-r", "--rxn_name")
    parser.add_argument("-q", "--quick", default=False, action="store_true")
    parser.add_argument("-n", "--num_processes")
    parser.add_argument("-f", "--force", default=False, action="store_true")
    parser.add_argument("-e", "--end_time")
    parser.add_argument("-p", "--prefix")

    args = parser.parse_args()

    if args.end_time is not None:
        print("Using custom end time of {}".format(args.end_time))
        end_t = float(args.end_time)
    if args.composition is not None:
        print("Using composition {}".format(args.composition))
        compositions = [args.composition]
    if args.rxn_name is not None:
        print("Using reaction {}".format(args.rxn_name))
        rxn_name = args.rxn_name
    if args.quick:
        print("Running in quick mode")
        Das = [d for d in Das if d <= 1e3]
        compositions = compositions if args.composition is None else [args.composition]
    if args.num_processes is not None:
        num_processes = int(args.num_processes)
    if args.force:
        save_output = True
        load_output = False
    if args.prefix is not None:
        prefix = args.prefix

##########################
# Ready output directory #
##########################

output_path = Path("figs",reference,rxn_name,prefix) if prefix is not None else Path("figs",reference,rxn_name)
output_path.mkdir(parents=True, exist_ok=True)
pickle_path = Path(output_path,"_outs.pickle")


######################################
# Plot pressure-temperature profiles #
######################################

fig = plt.figure(figsize=(5,7))
cmap1 = plt.cm.get_cmap("coolwarm_r")
ax1 = plt.gca()
ax1.set_prop_cycle(plt.cycler("color", cmap1(np.linspace(0., 1., len(tectonic_settings)))))

geotherm_latex_table = {
    "heading": ["" for s in tectonic_settings],
    "As": np.zeros(len(tectonic_settings)),
    "L0": np.zeros(len(tectonic_settings)),
    "z0": np.zeros(len(tectonic_settings)),
    "hr0": np.zeros(len(tectonic_settings)),
    "qs0": np.zeros(len(tectonic_settings)),
    "qs1": np.zeros(len(tectonic_settings)),
    "T0": np.zeros(len(tectonic_settings)),
    "T1": np.zeros(len(tectonic_settings)),
}

for num_setting, setting in enumerate(tectonic_settings):
    z0 = setting["z0"]
    L0 = setting["L0"]
    shortening = 1
    z1 = setting["z1"]
    hr0 = setting["hr0"]
    conductivity = setting["k"]
    Ts = setting["Ts"]
    Tlab = setting["Tlab"]
    As = setting["As"]

    depths = np.linspace(0,1)
    depths_sc = L0*depths
    P =  depths_sc * crustal_rho * gravity / 1e5
    T, qs0 = geotherm_steady(depths,
                        L0*shortening,
                        shortening,
                        Ts=Ts,
                        Tlab=Tlab,
                        k=conductivity,
                        A=As,
                        hr0=hr0)


    p = plt.plot(T-273.15, depths_sc/1e3,linewidth=1, alpha=0.5)
    color = plt.gca().lines[-1].get_color()
    #plt.plot(T[-1]-273.15, depths_sc[-1]/1e3, "x", color=color, )
    T0, _qs = geotherm_steady(z0/L0,
                    L0*shortening,
                    shortening,
                    Ts=Ts,
                    Tlab=Tlab,
                    k=conductivity,
                    A=As,
                    hr0=hr0)
    plt.plot(T0-273.15, z0/1e3,'.',color=color,alpha=1)
    
    shortening = z1/z0
    T, qs1 = geotherm_steady(depths,
                        L0*shortening,
                        shortening,
                        Ts=Ts,
                        Tlab=Tlab,
                        k=conductivity,
                        A=As,
                        hr0=hr0,)
    T1, _qs = geotherm_steady(z0/L0,
                L0*shortening,
                shortening,
                Ts=Ts,
                Tlab=Tlab,
                k=conductivity,
                A=As,
                hr0=hr0)
    plt.plot(T-273.15, depths_sc/1e3*shortening, "--", color=color,alpha=0.5,linewidth=1)

    Tts = np.zeros(100)
    zts = np.zeros(100)
    for i,t in enumerate(np.linspace(0,1,100)):
        s = 1 + (z1/z0 - 1)*t
        Tt, _q = geotherm_steady(z0/L0,
                L0*s,
                s,
                Ts=Ts,
                Tlab=Tlab,
                k=conductivity,
                A=As,
                hr0=hr0)
        Tts[i] = Tt
        zts[i] = z0*s
    
    label = setting["setting"].replace("hot","H").replace("transitional","T").replace("cold","C")
    
    plt.plot(np.array(Tts)-273.15,np.array(zts)/1.e3,'-',alpha=1,linewidth=1.2,color=color,label=label)

    # points after shortening
    plt.plot(T1-273.15, z0/1e3*shortening,'.',color=color)
    #plt.plot(T[-1]-273.15, depths_sc[-1]/1e3*shortening, "x", color=color)
    #plt.plot(T[0]-273.15, depths_sc[0]/1e3,'kx')
    geotherm_latex_table["heading"][num_setting] = setting["setting"].replace("hot","H").replace("cold","C").replace("transitional","T")
    geotherm_latex_table["As"][num_setting] = As*1.e6
    geotherm_latex_table["L0"][num_setting] = L0/1.e3
    geotherm_latex_table["z0"][num_setting] = z0/1.e3
    geotherm_latex_table["hr0"][num_setting] = hr0/1.e3
    geotherm_latex_table["qs0"][num_setting] = qs0*1.e3
    geotherm_latex_table["qs1"][num_setting] = qs1*1.e3
    geotherm_latex_table["T0"][num_setting] = T0-273.15
    geotherm_latex_table["T1"][num_setting] = T1-273.15

plt.legend()
ax1.set_ylabel("depth (km)")
ax1.set_xlabel("$T$ (°C)")

ax1.invert_yaxis()
plt.savefig(Path(output_path,"{}.{}".format("_geotherms", "pdf")), metadata=pdf_metadata)
plt.savefig(Path(output_path,"{}.{}".format("_geotherms", "png")))

ax1.set_ylim([0,120])
ax1.set_xlim([150,1100])
plt.savefig(Path(output_path,"{}.{}".format("_geotherms_inverted", "pdf")), metadata=pdf_metadata)
plt.savefig(Path(output_path,"{}.{}".format("_geotherms_inverted", "png")))

############################################
# Write LaTeX for 'tectonic setting' table #
############################################

table_body = """\\begin{{tabular}}{{cc{}}}
\\toprule
& & {} \\\\
\\midrule
$A_{{s}}$ & \\si{{\\uW\\per\\m\\cubed}} & {} \\\\
$L_{{0}}$ & \\si{{\\km}} & {} \\\\
$z_{{0}}$ & \\si{{\\km}} & {} \\\\
$h_{{r0}}$ & \\si{{\\km}} & {} \\\\
\\midrule
$q_{{s0}}$ & \\si{{\\mW\\per\\m\\squared}} & {} \\\\
$q_{{s1}}$ & \\si{{\\mW\\per\\m\\squared}} & {} \\\\
$T_{{0}}$ & \\si{{\\degreeCelsius}} & {} \\\\
$T_{{1}}$ & \\si{{\\degreeCelsius}} & {} \\\\
\\bottomrule
\\end{{tabular}}
""".format(
    "c"*len(geotherm_latex_table["heading"]),
    " & ".join(geotherm_latex_table["heading"]),
    " & ".join(["{:.2f}".format(val) for val in geotherm_latex_table["As"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["L0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["z0"]]),
    " & ".join(["{:.1f}".format(val) for val in geotherm_latex_table["hr0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["qs0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["qs1"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["T0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["T1"]]),

)

with open(Path(output_path,"_geotherms_table.tex"), "w") as fil:
    fil.writelines(table_body)

############################
# Setup initial conditions #
############################

# Create scenarios for all combinations 
# of setting, composition, and Da       
scenarios = []
for setting in tectonic_settings:
    for da in Das:
        for comp in compositions:
            scenario = setting.copy()
            scenario["Da"] = da
            scenario["composition"] = comp
            scenarios.append(scenario)

# Get equilibrium mantle densities
ipyrolite = get_rho_interpolator("xu_2008_pyrolite")
iharzburgite = get_rho_interpolator("xu_2008_harzburgite")

# Get reaction/phase/endmember parameters
rxn = get_reaction(rxn_name)
phase_names, endmember_names = get_names(rxn)
I = len(rxn.phases())
_Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(_Kis)

# Get equilibrated initial condition (with cache to avoid repeating work)
ic_cache = {}
def get_initial_composition(T0,P0,rxn,composition):
    slug = "{:.4f}-{:.4f}-{}-{}".format(T0,P0,rxn.name(),composition)
    if slug in ic_cache:
        return ic_cache[slug]

    # Get equilibrium composition at arbitrary P,T
    mi0, Xik0, phii0, Cik0 = get_point_composition(rxn, composition) 
    _Cik0 = x2c(rxn, Xik0) if Cik0 is None else Cik0
    _mi0 = phi2m(rxn, phii0, _Cik0) if mi0 is None else mi0

    # Equilibrate model at initial (T0, P0)
    ode = ScipyPDReactiveODE(rxn)
    ode.solve(T0, P0,_mi0,_Cik0,1,Da=1e5,eps=eps)
    
    rho0 = ode.final_rho()*100 # kg/m3
    Cik0 = ode.sol.y[ode.I:ode.I+ode.K,-1] # -1 = final time step
    mi0 = ode.sol.y[:ode.I,-1]

    ic_cache[slug] = (Cik0, mi0, rho0)
    return Cik0, mi0, rho0

# Get full initial conditions
def get_initial_state(rxn,scenario):

    composition = scenario['composition']
    L0 = scenario["L0"]
    z0 = scenario["z0"]
    As = scenario["As"]
    hr0 = scenario["hr0"]
    k = scenario["k"]
    Ts = scenario["Ts"]
    Tlab = scenario["Tlab"]
    z1 = scenario["z1"]

    # Initial temperature
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

    # Final temperature
    thickening = z1/z0
    T1, qs1 = geotherm_steady(
        z0/L0, # between 0 and 1
        L0*thickening,
        thickening=thickening, # gt 1
        Ts=Ts,
        Tlab=Tlab,
        k=k,
        A=As,
        hr0=hr0
    )

    # Initial pressure
    P0 = crustal_rho * gravity * z0/1e5 # bar

    # Initial composition (equilibrium at T0,P0)
    Cik0, mi0, rho0 = get_initial_composition(T0,P0,rxn,composition)

    return T0, P0, Cik0, mi0, rho0, qs0, T1, qs1

# Gets all initial conditions and save to scenario
def setup_ics(scenario):
    rxn = get_reaction(rxn_name)
    T0, P0, Cik0, mi0, rho0,qs0,T1,qs1 = get_initial_state(rxn,scenario)
    scenario["T0"] = T0
    scenario["T1"] = T1
    scenario["qs0"] = qs0
    scenario["qs1"] = qs1
    scenario["P0"] = P0
    scenario["Cik0"] = Cik0
    scenario["mi0"] = mi0
    scenario["rho0"] = rho0
    return scenario

#############################################
# Define RHS of differential system of eqns #
#############################################

def rhs(t,u,rxn,scale,Da,L0,z0,As,hr0,conductivity,T_surf,Tlab):

    # Damkholer number could be modulated here if wanted

    # Extract variables

    # NOTE: this is here called Fi for consistency with the paper. Elsewhere in the file it is called mi
    Fi = u[:I] # phase mass fractions
    
    # NOTE: this is here called cik for consistency with the paper. Elsewhere in the file it's called Cik
    cik = u[I:I+K] # endmember mass fractions

    T = u[I+K]
    P = u[I+K+1] 
 
    # Scalings
    Ts = scale["T"]*T
    Ps = scale["P"]*P
    rho0 = scale["rho"]
    dz = scale["h"] # length scale h0

    # reshape C
    C = rxn.zero_C() # object with correct shape
    Kis = 0
    for i,Ki in enumerate(_Kis):
        C[i] = cik[Kis:Kis+Ki]
        Kis = Kis+Ki

    # regularize C by taking max of C and eps
    Cs = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    Cs = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    rhoi = np.array(rxn.rho(Ts, Ps, Cs)) # phase densities $\rho_i$
    V = np.sum(Fi/rhoi) # total volume
    
    # regularize F by adding eps
    Fis = np.asarray(Fi)
    Fis = Fi + eps

    # Get dimensionless Gammas from reaction
    Gammai = np.asarray(rxn.Gamma_i(Ts,Ps,Cs,Fi))
    gamma_ik = rxn.Gamma_ik(Ts,Ps,Cs,Fi)
    Gammaik = np.zeros(K)
    sKi = 0
    for i in range(I):
        for k in range(_Kis[i]):
            Gammaik[sKi+k] = gamma_ik[i][k]
        sKi += _Kis[i]
    
    # Calculate reactive mass flux (scaled Gammas)
    du = np.zeros(u.shape)
    sKi = 0
    for i in range(I):
        du[i] = Da*rho0*Gammai[i]*V
        for k in range(_Kis[i]):
            GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
            du[I+sKi+k] = Da*rho0*GikcGi*V/Fis[i]
        sKi += _Kis[i]

    # dT/dt
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

    # dP/dt
    dP = crustal_rho * gravity / 1e5 * dz # bar 
    dP = dP/scale["P"]

    # add dT and dP to outputs
    du[I+K:] = np.array([dT, dP])

    return du

############################################
# Function that runs scenarios in parallel #
############################################

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

    # Set reaction's characteristic Arrhenius temperature (T_r in paper)
    rxn.set_parameter("T0",2000+273.15) # 2000 C gives a good fit to data

    scale= {"T":T0, "P":P0, "rho":rho0, "h":(z1-z0)}

    # Set up vector of initial conditions
    u0=np.empty(I+K+2) # each endmember, each phase, pressure, temperature
    u0[:I] = mi0 # intial phase mass fractions
    u0[I:I+K] = Cik0 # initial endmember mass fractions
    u0[I+K:] = np.array([1., 1.]) # scaled T0 and P0 (T/T0 = P/P0 = 1)

    #event = lambda t,y,rxn,scale,Da: Pmax - y[-1]*scale["P"] # when pressure = Pmax
    #event.terminal = True

    # Rescale damkohler number in case the end time is not 1
    _da = Da * end_t
    args = (rxn,scale,_da,L0,z0,As,hr0,k,Ts,Tlab)

    # Solve IVP using BDF method
    sol = solve_ivp(rhs, [0, end_t], u0, args=args, dense_output=True, method="BDF", rtol=rtol, atol=atol, events=None)
    
    # resample solution
    t = np.linspace(0,end_t,1000)
    y = sol.sol(t)

    # Dimensionalize back T,P
    T = y[-2]*scale["T"] # K
    P = y[-1]*scale["P"] # bar

    mi_times  = y[:I].T # vector for each timestep
    Cik_times = y[I:I+K].T # 2d array for each timestep

    # reshape C
    def reshape_C(Cik):
        C = rxn.zero_C()
        k = 0
        for i,Ki in enumerate(_Kis):
            C[i] = Cik[k:k+Ki]
            k = k+Ki
        return C
    
    Cs_times = [reshape_C(Cik) for Cik in Cik_times] # vector for each timestep

    # Calculate rho for each timestep as 1/sum_i(F_i/rho_i)
    # for which we need the endmember compositions as a vector for each phase (Cs_times)
    rho = [1/sum(mi_times[t]/rxn.rho(T[t], P[t], Cs))/10 for t,Cs in enumerate(Cs_times)]

    print("{} P_end = {:.2f} Gpa. T_end = {:.2f} K. Used {:n} steps: ".format(sol.message,P[-1]/1e4,T[-1],len(sol.t)))

    # Back-calculate depth from pressure
    depth_m = (P*1e5) / gravity / crustal_rho
    
    scenario["T"] = T # K
    scenario["P"] = P # bar
    scenario["rho"] = rho # g/cm3
    scenario["mi"] = mi_times # phase mass fractions
    scenario["Cik"] = Cik_times # endmember mass fractions
    scenario["Xik"] = np.asarray([rxn.C_to_X(c) for c in Cs_times], dtype="object") # endmember mol. fractions
    scenario["z"] = depth_m
    scenario["time"] = t
    return scenario

#####################################
# Deal with result loading & saving #
#####################################

if load_output:
    print("Loading pickle from file file {}".format(pickle_path))
    with open(pickle_path, 'rb') as pickle_file:
        scenarios_out = pickle.load(pickle_file)
else:
    print("Preparing to run {} scenarios".format(len(scenarios)))
    scenarios_in = [setup_ics(s) for s in scenarios]

    # run for varying damkhoeler numbers
    with Pool(num_processes) as pool:
        # blocks until all finished
        scenarios_out = pool.map(run_experiment, scenarios_in)

if save_output:
    with open(pickle_path, 'wb') as pickle_file:
        pickle.dump(scenarios_out, pickle_file)

###################
# Post processing #
###################

# this needs to be in global scope
bin_widths = {"1Myr":1*Myr, "100kyr":100*kyr, "50kyr":50*kyr, "10kyr":10*kyr}

for out in scenarios_out:
    # Process each scenario output
    rho = np.array(out["rho"]) # g/cm3
    depth_m = out["z"]
    T = out["T"] # K
    t = out["time"] # unitless
    P = out["P"] # bar
    z1 = out["z1"]
    z0 = out["z0"]
    
    rho_pyrolite = ipyrolite((T, P/1e4))/10

    max_rho = np.nanmax(rho)
    max_rho_py = np.nanmax(rho_pyrolite)

    # find the critical depth (pressure, temperature, etc)
    if max_rho > max_rho_py:
        # find the greatest depth at which rho exceeds pyrolite
        critical_indices = [i for i,r in enumerate(rho) if r > rho_pyrolite[i]]
        if len(critical_indices) == 0:
            # never goes critical
            critical_depth = 85.e3 # approx. coesite transition?, can change this later
            critical_pressure = 30e3 # bar
            critical_temperature = T[-1] # K
            critical_time = end_t+1
        elif len(critical_indices) == len(rho):
            # always critical
            critical_depth = depth_m[0]
            critical_pressure = P[0] # bar
            critical_temperature = T[0] # K
            critical_time = t[0]
        else:
            # interesting 
            first_critical_index = critical_indices[0]
            critical_index = first_critical_index
            critical_depth = depth_m[critical_index]
            critical_pressure = P[critical_index] # bar
            critical_temperature = T[critical_index] # K
            critical_time = t[critical_index]
    else:
        critical_depth = 85.e3 # approx. coesite transition?, can change this later
        critical_pressure = 30e3 # bar
        critical_temperature = T[-1] # K
        critical_time = end_t+1
    out["critical_depth"] = critical_depth
    out["critical_pressure"] = critical_pressure
    out["critical_temperature"] = critical_temperature
    out["critical_time"] = critical_time

    # Average density contrast of any unstable root
    delta_rho = rho - rho_pyrolite
    delta_rho[delta_rho < 0] = np.nan
    effective_delta_rho = np.nanmean(delta_rho)
    out["effective_delta_rho"] = effective_delta_rho
    out["max_rho"] = max_rho

    # Densification rate
    time = np.linspace(0,end_t, rho.size) * t0 # seconds
    time_Myr = time/Myr
    densification_rate = np.diff(rho*1000)/np.diff(time_Myr) # kg/m3/Myr
    densification_rate = np.insert(densification_rate, 0, 0.)
    out["densification_rate"] = densification_rate
    out["time_Myr"] = time_Myr
    for i,phase in enumerate(phase_names):
        phase_mis = out["mi"][:,i] # wt%
        if(phase=='Feldspar'):
            plag_frac = phase_mis.copy()
            # find the index at which plagioclase drops below 1 wt%
            plag_out_indices = [i for i,X in enumerate(plag_frac) if X < 0.025]
            if len(plag_out_indices) == 0:
                out['plag_out_depth'] = np.nan 
                out['plag_out_pressure'] = np.nan # bar
                out['plag_out_temperature'] = np.nan # K
            else:
                first_plag_out = plag_out_indices[0]
                plag_out_index = first_plag_out
                out['plag_out_depth'] = depth_m[plag_out_index]
                out['plag_out_pressure'] = P[plag_out_index] # bar
                out['plag_out_temperature'] = T[plag_out_index] # K

    # Max densification rate, binned and averaged over time bins,
    # limited to the P-T space of the plag-out reaction (not the garnet-in reaction)
    for bin_width_string, bin_width in bin_widths.items():
        bins = np.arange(0., t0, int(bin_width))
        digitized_t = np.digitize(time, bins) # assigns the index of a bin to each point
        plag_out_mask = P/1.e4 < 0.5 + 1./1000.*(T-273.15-300.) # NOTE: 'True' signifies masking OUT (nan)
        dens_rate_masked = ma.masked_array(densification_rate,mask=plag_out_mask)
        bin_means = [dens_rate_masked[digitized_t == i].mean() for i in range(len(bins))]
        out["max_densification_rate_"+bin_width_string] = np.nanmax(bin_means)

# Create '_critical' plot
selected_compositions = [
    #"mackwell_1998_maryland_diabase_norm",
    #"hacker_2015_md_xenolith_norm",
    #"sammon_2021_lower_crust_norm",
    #"sammon_2021_deep_crust_norm",
    "si50_norm",
    "si51_norm",
    "si52_norm",
    "si53_norm",
    "si54_norm",
    "si55_norm",
    "si56_norm",
    "si57_norm",
    "si58_norm"
]

selected_outputs = [o for o in scenarios_out if o["composition"] in selected_compositions]

depths = [o['critical_depth'] for o in selected_outputs]

T0s = np.asarray([o["T0"] - 273.15 for o in selected_outputs])

def T2size(T):
    return ((1 - (T-np.min(T))/(np.max(T)-np.min(T)))*6 + 3)**2
def size2T(s):
    return -((np.sqrt(s)-3)/6 - 1)*(np.max(T0s)-np.min(T0s))+np.min(T0s)

sizes = T2size(T0s)

comps = [o["composition"] for o in selected_outputs]
codes = [selected_compositions.index(c) for c in comps]

fig = plt.figure(figsize=(5,7))
s = plt.scatter([o["Da"] for o in selected_outputs], np.asarray(depths)/1e3, s=sizes, c=codes, cmap="tab10", alpha=0.5)

# produce a legend with the unique colors from the scatter
ax = plt.gca()
handles, labels = s.legend_elements(
    prop="colors",
    #func=lambda x: comps[int(x)].capitalize().replace("_", " ")
)
legend1 = ax.legend(handles,
                    [c.capitalize().replace("_"," ") for c in selected_compositions],
                    loc="lower left",
                    bbox_to_anchor=(1.04,0),
                    title="$X_{{bulk}}$")

ax.add_artist(legend1)
# produce a legend with a cross-section of sizes from the scatter
handles, labels = s.legend_elements(
    prop="sizes",
    alpha=0.6,
    fmt="{x:.0f}",
    func=size2T
)

legend2 = ax.legend(handles, 
                    labels, 
                    loc="upper left",
                    bbox_to_anchor=(1.04,1),
                    title="$T_0$ (°C)")
ax.set_ylim(top=30, bottom=85)
plt.semilogx()
plt.xlabel("Da")
plt.ylabel("Critical depth (km)")
plt.savefig(Path(output_path,"{}.{}".format("_critical", "pdf")), metadata=pdf_metadata, bbox_extra_artists=(legend1,legend2), bbox_inches='tight')
plt.savefig(Path(output_path,"{}.{}".format("_critical", "png")),bbox_extra_artists=(legend1,legend2), bbox_inches='tight')

# Write summary data to CSV
with open(Path(output_path,'_critical.csv'),'w') as csvfile:
    fieldnames = [
        'setting',
        'L0',
        'z0',
        'z1',
        'As',
        'hr0',
        'k',
        'Ts',
        'Tlab',
        'Da',
        'composition',
        'T0',
        'T1',
        'P0',
        'rho0',
        'critical_depth',
        'critical_pressure',
        'critical_temperature',
        'plag_out_pressure',
        'plag_out_temperature',
        'plag_out_depth',
        'qs0',
        'qs1',
        'effective_delta_rho',
        'max_rho'
    ]
    for key,val in bin_widths.items():
        fieldnames.append('max_densification_rate_'+key)
    writer = csv.DictWriter(csvfile, fieldnames,extrasaction='ignore')
    writer.writeheader()
    for out in scenarios_out:
        writer.writerow(out)
    print("wrote _critical.csv")

# Output plots for every composition & tectonic setting
for composition in compositions:
    for tectonic_setting in tectonic_settings:
        setting = tectonic_setting["setting"]

        # Find outputs with same setting & same composition,
        # and sort by Da
        outs_c = sorted([out for out in scenarios_out if (out["composition"] == composition and out["setting"] == setting)], key=lambda out: out["Da"])
        
        # grab 'constants' from first output
        base = outs_c[0] 
        T = base["T"]
        P = base["P"]
        rho0 = base["rho0"]
        L0 = base["L0"]
        z0 = base["z0"]
        z1 = base["z1"]

        _Das = [o["Da"] for o in outs_c]
        
        # Back-calculate reaction coefficient from Da
        r0 = [da*rho0/t0 for da in _Das] # Gamma0 (kg/m3/s)
        S0 = 6000 # 1/m
        reaction_rate_per_surface = [r/S0 for r in r0] # r0 (kg/m2/s)
        reaction_rate_per_surface_gcm = [r*1000/100/100*yr for r in reaction_rate_per_surface] # g/cm2/yr

        # Setup subplots and figure
        num_subplots = 3 + len(phase_names) + 2
        subplot_mosaic = [["rho","rho"] + phase_names + ["An", "Jd", "T"]] # axes are named
        fig = plt.figure(figsize=(3*num_subplots,12))
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        axes = fig.subplot_mosaic(subplot_mosaic)

        # Invert y axis because it represents depth
        [ax.invert_yaxis() for label,ax in axes.items()]

        # Setup color cycling through all Das
        num_lines = len(_Das)
        greys = plt.cm.get_cmap("Greys")
        axes["rho"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))
        axes["An"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))
        axes["Jd"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))

        # Plot density vs depth
        for i, obj in enumerate(outs_c):
            ax = axes["rho"]
            ax.plot(obj["rho"], obj["z"])

        rho_pyrolite=ipyrolite((T, P/1e4))
        rho_harzburgite=iharzburgite((T,P/1e4))

        ax.plot(rho_pyrolite/10, outs_c[0]["z"], "r:")
        ax.plot(rho_harzburgite/10, outs_c[0]["z"], "g:")
        ax.legend(["$Da = ${:.1e}".format(d) for d in _Das] +  (["pyrolite", "harzburgite"]), loc="upper right")

        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Density")
        ax.set_xlim([2.8, 3.8])

        # Plot temperature vs depth
        ax = axes["T"]
        ax.plot(T-273.15, base["z"]/1e3, linewidth=2)
        ax.set_xlabel("T (°C)")
        ax2 = ax.twinx()
        ax2.plot(T-273.15, P/1e4, alpha=0)
        ax2.invert_yaxis()
        ax2.set_ylabel("P (GPa)")

        # Plot all phase mis (mass fractions)
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

        # get Plagioclase anorthite composition
        ax = axes["An"]
        ax.set_xlim([0., 1.0])
        ax.set_ylabel(None)
        ax.set_yticklabels([])
        ax.set_xlabel("$X_{\mathrm{An}}$")

        # get Cpx jadeite content
        ax = axes["Jd"]
        ax.set_xlim([0., 1.0])
        ax.set_ylabel(None)
        ax.set_yticklabels([])
        ax.set_xlabel("$X_{\mathrm{Jd}}$")

        # Plot X_An and X_Jd vs depth
        for j, obj in enumerate(outs_c):
            XAn = [x[3][0] for x in obj["Xik"]]
            XJd = [x[0][4] for x in obj["Xik"]]

            axes["An"].plot(XAn, obj["z"])
            axes["Jd"].plot(XJd, obj["z"])

        fig.suptitle("{}, {}, $v_0=${:.1f} km/Myr, $S0=${} 1/m".format(composition.capitalize().replace("_"," "), setting.replace("_",", "), moho_descent_rate/1e3*yr*1e6, S0),y=0.9)
        plt.savefig(Path(output_path,"{}.{}.{}".format(setting,composition,"pdf")), metadata=pdf_metadata)
        plt.savefig(Path(output_path,"{}.{}.{}".format(setting,composition,"png")))
        plt.close(fig)

for tectonic_setting in tectonic_settings:

        setting = tectonic_setting["setting"]
        outs_c = sorted([out for out in scenarios_out if (out["composition"] in selected_compositions and out["setting"] == setting)], key=lambda out: out["Da"])
        
        # grab 'constants' from first output
        base = outs_c[0] 
        T = base["T"]
        P = base["P"]
        # Setup figure for rho-profile mosaic
        fig = plt.figure(figsize=(len(selected_compositions)*2,6.25))
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        axes = fig.subplot_mosaic([selected_compositions])
        [ax.set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines)))) for label,ax in axes.items()]

        # Invert y axis because it represents depth
        [ax.invert_yaxis() for label,ax in axes.items()]
        [ax.set_xlim([2.8,3.5]) for label,ax in axes.items()]
        [ax.set_xticks([2.8,3.0,3.2,3.4]) for label,ax in axes.items()]
        [ax.set_ylim([80,30]) for label,ax in axes.items()]
        rho_pyrolite=ipyrolite((T, P/1e4))/10
        rho_harzburgite=iharzburgite((T,P/1e4))/10
        for i, obj in enumerate(outs_c):
            ax = axes[obj["composition"]]
            ax.plot(obj["rho"], obj["z"]/1e3)
            if obj["Da"] == 1:
                ax.plot(rho_pyrolite, obj["z"]/1e3, "r:")
                ax.plot(rho_harzburgite, obj["z"]/1e3, "g:")
        plt.savefig(Path(output_path,"_collage.{}.{}".format(setting,"pdf")), metadata=pdf_metadata)
        plt.savefig(Path(output_path,"_collage.{}.{}".format(setting,"png")))
        plt.close(fig)
