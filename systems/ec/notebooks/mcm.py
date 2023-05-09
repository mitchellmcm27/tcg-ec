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

        self.Cik0 = Cik0 if Cik0 is not None else self.x2c(self.rxn, Xik0)
        self.mi0 = mi0 if mi0 is not None else self.phi2m(
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

    def x2c(self, rxn, Xik0):
        return np.asarray([c for (i, ph) in enumerate(rxn.phases()) for c in ph.x_to_c(Xik0[i])])

    def phi2m(self, rxn, vi0, Cik0, T=None, P=None):
        '''Converts phase modes in volume fraction to mass fraction given an intial EM composition in mass fractions.'''
  
        densities = []
        i = 0
        if T is None:
            T = 300 # K
        if P is None:
            P = 0.0001013 # GPa

        for ph in rxn.phases():
            n = len(ph.endmembers())
            densities.append(ph.rho(T, GPa2Bar(P), Cik0[i:i+n]))
            i = i+n
        #print(densities)
        mass = np.sum(np.asarray(densities) * np.asarray(vi0))
        #print(mass)
        
        mi0 = np.asarray([v*densities[i]/mass for (i, v) in enumerate(vi0)])
        print(np.sum(mi0))
        return mi0

    @staticmethod
    def get_reaction(rxnName):
        pv = repr(sys.version_info.major)+'.'+repr(sys.version_info.minor)
        path = os.path.join(os.path.pardir, 'database', 'install', rxnName,
                            'lib', 'python'+pv, 'site-packages/')  # the final slash is necessary!
        sys.path.append(path)
        #print(path)
        tcgdb = __import__('py_'+rxnName)
        func = getattr(tcgdb, rxnName)
        rxn = func()  # <-- this should work!
        return rxn

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
        
        bdfdiag = decorate(PDReactiveProfileDiagnostics)(self.rxn, bdfgrid)

        s2 = bdfdiag.plot_rho_contours()

        if plot_phases:
            bdfdiag.plot_phases()
        bdfdiag.plot_path()
        return bdfdiag