import sys, os
sys.path.append(os.path.join(os.path.pardir, 'python'))

import math
import pandas as pd
import numpy as np
from pathlib import Path
import pickle
from tcg_slb.base import *
from tcg_slb.phasediagram.base import PDReactiveGrid, PDReactiveProfile, PDReactiveGridDiagnostics, PDReactiveProfileDiagnostics
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
import sympy as sym
sym.init_printing()

class EcModel():
    Tmin = 0.
    Tmax = 0.
    nT = 0
    T = None

    Pmin = 0.
    Pmax = 0.
    nP = 0
    P = None

    rxn = None

    Cik0 = None
    mi0 = None

    # results
    grid = None
    bdfdiag = None
    ode = None

    def __init__(self, referenceName, rxnName, Tmin=773., Tmax=1273., nT=60, nP=60, Pmin=0.5, Pmax=2.5, Cik0=None, mi0=None, phii0=None, Xik0=None, T0=None, P0=None, domain="grid", **kwargs):
        self.referenceName = referenceName
        self.rxnName = rxnName

        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nT = nT
        self.nP = nP
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.T0 = T0
        self.P0 = P0
        self.domain = domain
        self.rxn = self.get_reaction(rxnName)
        self.phaseNames = [ph.name() for ph in self.rxn.phases()]

        if self.domain=="grid" and nT > 0 and nP > 0:
            self.T = np.linspace(Tmin, Tmax, nT)
            self.P = np.linspace(Pmin, Pmax, nP)
        if self.domain=="profile" and nT> 0:
            num_points = math.floor(nT)
            self.T = np.linspace(Tmin,Tmax,num_points)
            self.P = np.linspace(Pmin,Pmax,num_points)

        self.Cik0 = Cik0 if Cik0 is not None else x2c(self.rxn, Xik0)
        self.mi0 = mi0 if mi0 is not None else phi2m(
            self.rxn, phii0, self.Cik0, P=P0, T=T0)

        #print(self.Cik0)
        #print('')
        #if phii0 is not None:
        #    print(phii0)
        #    print(sum(phii0))
        #    print('')
        #print(self.mi0)
        #print(sum(self.mi0))

    def run(self, end_t=1e2, reload=False, save=False, plot=True, plot_phases=True, **kwargs):
        # initial temperature, pressure and phase volume fraction
        Ti = self.T0 # if self.T0 else 900. # kelvin
        pi = self. P0 # GPa2Bar(self.P0 if self.P0 else 1.25)  # bars
        if(Ti is not None and pi is not None):
            Ti = self.T0
            pi = GPa2Bar(self.P0)
            ode = ScipyPDReactiveODE(self.rxn)
            self.ode = ode
            ode.solve(Ti, pi, self.mi0, self.Cik0, end_t, **kwargs)
            #print(ode.stime)
            #print(ode.final_phases(1.e-2))
            if plot:
                ode.plot()
                self.display_initial_final()

        if self.domain == "grid" and self.T is not None and self.P is not None:
            grid = self.solve_reaction_grid(
                reload=reload, save=save, end_t=end_t, **kwargs)
            self.grid = grid
            if(plot):
                bdfdiag = self.plot_reaction_grid(grid, plot_phases=plot_phases)
                self.bdfdiag = bdfdiag
            else:
                bdfdiag = None
            return self.rxn, grid, ode, bdfdiag
        elif self.domain=="profile" and self.T is not None and self.P is not None:
            grid = self.solve_reaction_profile(reload=reload, save=save, end_t=end_t, **kwargs)
            print("done solving profile")
            self.grid = grid
            if(plot):
                bdfdiag = self.plot_reaction_profile(grid, plot_phases=plot_phases)
                self.bdfdiag = bdfdiag
            else:
                bdfdiag = None
            return self.rxn, grid, ode, bdfdiag
        else:
            return self.rxn, None, ode, None
        
    def mi_final(self):
        return self.ode.sol.y[:self.ode.I, -1]
    def Cik_final(self):
        Cik = self.ode.sol.y[self.ode.I:self.ode.I+self.ode.K].T
        Cik_end = Cik[-1]
        Cs = self.ode.regularizeC(Cik_end)
        return [cik for ci in Cs for cik in ci]
    
    def display_initial_final(self):
        mi = self.mi_final()
        table = [self.ode.mi0, mi]
        # convert to volume!
        phaseNames = [ph.name() for ph in self.rxn.phases()]
        df = pd.DataFrame(
            table, 
            columns=phaseNames, 
            index=['Wt% (initial)', 'Wt% (final)']
        )

        try:
            from IPython.display import display, HTML
            display(df)
        except:
            print(df)

    def get_pickle_path(self):
        if(self.domain=="profile"):
            return Path('output', self.referenceName, self.rxn.name() + '_profile.pickle')
        return Path('output', self.referenceName, self.rxn.name() + '.pickle')

    def load_grid(self):
        filename = self.get_pickle_path()
        with open(filename, 'rb') as pfile:
            bdfgrid = pickle.load(pfile)
        self.bdfgrid = bdfgrid
        return bdfgrid

    def save_grid(self, bdfgrid):
        filename = self.get_pickle_path()
        filename.parent.mkdir(exist_ok=True, parents=True)
        with open(filename, 'wb+') as pfile:
            pickle.dump(bdfgrid, pfile)

    def solve_reaction_profile(self, end_t=1e2, reload=False, save=False, odeClass=ScipyPDReactiveODE, **kwargs):
        if reload:
            return self.load_grid()
        i0 = 2  # doesn't matter as long as you pass cik0
        bdfgrid = PDReactiveProfile()
        bdfgrid.solve(self.rxn, odeClass, i0, [
                      'T', 'p'], self.T, self.P, end_t, Cik0=self.Cik0, mi0=self.mi0)
        if save:
            self.save_grid(bdfgrid)
        return bdfgrid
    
    def solve_reaction_grid(self, end_t=1e2, reload=False, save=False, odeClass=ScipyPDReactiveODE):
        if reload:
            return self.load_grid()

        i0 = 2  # doesn't matter as long as you pass cik0

        bdfgrid = PDReactiveGrid()

        bdfgrid.solve(self.rxn, odeClass, i0, [
                      'T', 'p'], self.T, self.P, end_t, Cik0=self.Cik0, mi0=self.mi0)

        if save:
            self.save_grid(bdfgrid)

        return bdfgrid

    def plot_reaction_grid(self, bdfgrid, plot_phases=True, figure_background=None):
        import matplotlib.pyplot as plt

        figure_xlim = [self.Tmin, self.Tmax]
        figure_ylim = [self.Pmin, self.Pmax]

        def decorate(pdrgd):
            def new_setup_axes(self, axi):
                if (figure_background is not None):
                    img = plt.imread(figure_background)
                    ip = axi.imshow(img)
                axi.axis('off')
                ax = axi.inset_axes([0.001, 0.006, 0.995, 0.991])
                ax.patch.set_alpha(0.0)

                ax.set_xlabel("Temperature (K)")
                ax.set_ylabel("Pressure (GPa)")
                ax.set_xlim(figure_xlim)
                ax.set_ylim(figure_ylim)
                return ax

            # replace the display with newdisplay
            pdrgd.setup_axes = new_setup_axes

            # return the modified student
            return pdrgd

        bdfdiag = decorate(PDReactiveGridDiagnostics)(self.rxn, bdfgrid)
        # s=bdfdiag.plot_rho()
        # s.set_clim([25., 38.])
        # s.set_cmap('jet')

        s2 = bdfdiag.plot_rho_contours()

        if plot_phases:
            bdfdiag.plot_phases()
        return bdfdiag

    def plot_reaction_profile(self, bdfgrid, plot_phases=True, figure_background=None):
        print("plot profile")
        import matplotlib.pyplot as plt

        figure_xlim = [self.T[0], self.T[-1]]
        def decorate(pdrgd):
                def new_setup_axes(self, ax):
                    ax.set_xlabel("Temperature (K)")
                    ax.set_xlim(figure_xlim)   
                    return ax

                # replace the display with newdisplay
                pdrgd.setup_axes = new_setup_axes

                # return the modified student
                return pdrgd
        
        print("Decorate")
        bdfdiag = decorate(PDReactiveProfileDiagnostics)(self.rxn, bdfgrid)

        print("plot_rho_contours")
        s2 = bdfdiag.plot_rho_contours()

        if plot_phases:
            print("plot_phases")
            bdfdiag.plot_phases()
        print("plot_path")
        bdfdiag.plot_path()
        return bdfdiag

def x2c(rxn, Xik0):
    return np.asarray([c for (i, ph) in enumerate(rxn.phases()) for c in ph.x_to_c(Xik0[i])])

def phi2m(rxn, phii0, Cik0, T=900.,p=10000.,eps=1e-5):
    densities = []
    C = rxn.zero_C()
    Ki = 0
    for i,ph in enumerate(rxn.phases()):
        n = len(ph.endmembers())
        C[i] = Cik0[Ki:Ki+n]
        Ki = Ki+n
    # regularize C
    C = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    C = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    densities = [ph.rho(T, p, C[i]) for i,ph in enumerate(rxn.phases())]
    mass_tot = np.sum(np.asarray(densities) * np.asarray(phii0))
    mi0 = np.asarray([v*densities[i]/mass_tot for (i, v) in enumerate(phii0)])
    return mi0

def get_reaction(rxnName):
    pv = repr(sys.version_info.major)+'.'+repr(sys.version_info.minor)
    path = os.path.join(os.path.pardir, 'database', 'install', rxnName,
                        'lib', 'python'+pv, 'site-packages/')  # the final slash is necessary
    sys.path.append(path)
    #print(path)
    tcgdb = __import__('py_'+rxnName)
    importer = getattr(tcgdb, rxnName)
    rxn = importer()
    return rxn

def get_names(rxn):
    phase_names = [p.name() for p in rxn.phases()]
    phase_names = [s.replace("_slb_ph","") for s in phase_names]

    endmember_names = [em.name() for p in rxn.phases() for em in p.endmembers()]
    endmember_names = [s.replace("_slb_em","") for s in endmember_names]
    return phase_names, endmember_names

def latex_reactions(rxn):
    '''
    rudnick_2014_lower_crust
    Reaction object: eclogitization_agu18_stx21_rx

    Phase 0 Clinopyroxene_slb_ph (cpx)
         Endmember 0 Diopside_slb_em : CaMgSi2O6_(cpx)
         Endmember 1 Hedenbergite_slb_em : CaFeSi2O6_(cpx)
         Endmember 2 Clinoenstatite_slb_em : Mg2Si2O6_(cpx)
         Endmember 3 CaTschermaks_slb_em : CaAl2SiO6_(cpx)
         Endmember 4 Jadeite_slb_em : NaAlSi2O6_(cpx)
    Phase 1 Orthopyroxene_slb_ph (opx)
         Endmember 0 Enstatite_slb_em : Mg2Si2O6_(opx)
         Endmember 1 Ferrosilite_slb_em : Fe2Si2O6_(opx)
         Endmember 2 MgTschermaks_slb_em : MgAl2SiO6_(opx)
         Endmember 3 OrthoDiopside_slb_em : CaMgSi2O6_(opx)
    Phase 2 Quartz_slb_ph (qtz)
         Endmember 0 Quartz_slb_em : SiO2_(qtz)
    Phase 3 Feldspar_slb_ph (plg)
         Endmember 0 Anorthite_slb_em : CaAl2Si2O8_(plg)
         Endmember 1 Albite_slb_em : NaAlSi3O8_(plg)
    Phase 4 Garnet_slb_ph (gt)
         Endmember 0 Pyrope_slb_em : Mg3Al2Si3O12_(gt)
         Endmember 1 Almandine_slb_em : Fe3Al2Si3O12_(gt)
         Endmember 2 Grossular_slb_em : Ca3Al2Si3O12_(gt)
         Endmember 3 MgMajorite_slb_em : Mg4Si4O12_(gt)
         Endmember 4 NaMajorite_slb_em : Na2Al2Si4O12_(gt)
    Phase 5 Kyanite_slb_ph (ky)
         Endmember 0 Kyanite_slb_em : Al2SiO5_(ky)

    Reaction 0
         0.666667 CaFeSi2O6_(cpx) + 0.333333 Mg2Si2O6_(opx) -> 0.666667 CaMgSi2O6_(cpx) + 0.333333 Fe2Si2O6_(opx)
    Reaction 1
         0.6 Fe2Si2O6_(opx) + 0.4 Mg3Al2Si3O12_(gt) -> 0.6 Mg2Si2O6_(opx) + 0.4 Fe3Al2Si3O12_(gt)
    Reaction 2
         0.75 CaFeSi2O6_(cpx) + 0.25 Mg3Al2Si3O12_(gt) -> 0.75 CaMgSi2O6_(cpx) + 0.25 Fe3Al2Si3O12_(gt)
    Reaction 3
         Mg2Si2O6_(opx) -> Mg2Si2O6_(cpx)
    '''
    
    import io
    from contextlib import redirect_stdout

    with io.StringIO() as buf, redirect_stdout(buf):
        rxn.report()
        report = buf.getvalue()

    lines = report.splitlines()
    ref = lines[0]
    name = lines[1]
    phases=[]
    reactions=[]
    names = []
    for i in range(2,len(lines)):
        line = lines[i]
        if(line[0:5]=="Phase"):
            phases.append(line[8:].strip())
        if(line[0:8]=="Reaction"):
            n = int(line[8:].strip())+1
            names.append("\\textbf{"+str(n)+".}")
            rxn = lines[i+1].strip()
            rxn = rxn.replace("0.333333","1/3")
            rxn = rxn.replace("0.666667","2/3")
            rxn = rxn.replace("0.75","3/4")
            rxn = rxn.replace("0.5","1/2")
            rxn = rxn.replace("0.25","1/4")
            rxn = rxn.replace("0.2","1/5")
            rxn = rxn.replace("0.4","2/5")
            rxn = rxn.replace("0.6","3/5")
            rxn = rxn.replace("0.8","4/5")
            rxn = rxn.replace("0.166667","1/6")
            rxn = rxn.replace("_","")
            rxn = rxn.replace("->", "=")
            rxn = "\ce{" + rxn + "}"
            reactions.append(rxn)

    # 1 & $\ce{2/3 CaFeSi2O6}^\text{cpx} + \ce{1/3 Mg2Si2O6}^\text{opx} = \ce{2/3 CaMgSi2O6}^\text{cpx} + \ce{1/3 Fe2Si2O6}^\text{opx}$ \\
    table = "\n".join(["{} & {} \\\\".format(names[i], reactions[i]) for i in range(len(names))])
    return table

def composition_to_label(c):
    default = c.replace("_", " ").capitalize()
    dic = {
        "hacker_2015_md_xenolith": "Hacker et al. (2015) median xenolith",
        "hacker_2015_bin_1": "Hacker et al. (2015) bin-1",
        "hacker_2015_bin_2": "Hacker et al. (2015) bin-2",
        "hacker_2015_bin_3": "Hacker et al. (2015) bin-3",
        "hacker_2015_bin_4": "Hacker et al. (2015) bin-4",
        "sammon_2021_lower_crust": "Sammon & McDonough (2021) lower crust",
        "sammon_2021_deep_crust": "Sammon & McDonough (2021) deep crust",
        "xu_2008_basalt": "Xu et al. (2008) basalt",
        "xu_2008_pyrolite": "Xu et al. (2008) pyrolite",
        "xu_2008_harzburgite": "Xu et al. (2008) harzburgite",
        "zhang_2022_cd07-2": "Zhang et al. (2022) granulite",
        "zhang_2006_mafic_granulite": "Zhang et al. (2006) granulite",
        "bhowany_2018_hol2a": "Bhowany et al. (2018) eclogite (hol2a)"
    }
    return dic.get(c, default)
