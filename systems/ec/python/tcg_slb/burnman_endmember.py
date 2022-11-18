import sympy as sym
import re
import types
from thermocodegen.coder import coder
from .endmember import TCGEndmember

class BurnmanEndmember(TCGEndmember):
    name    = None
    formula = None
    model   = None
    param_strs = None
    param_unit = None
    param_vals = None
    A = None # Helmholtz
    G = None # Gibbs
    syms = None
    param_syms = None
    model_type = 'TP' # or 'TV for Helmholtz/Stixrude'

    def __init__(self,name,formula,reference,**kwargs):
        super().__init__(name,formula,reference,'TP')

        self.param_strs = [ 'H_0', 'S_0','P_0','a_0','T_einstein','n','T_0',
                            'Cp_0', 'Cp_1', 'Cp_2', 'Cp_3',
                            'V_0','a_0',
                            'K_0','Kprime_0','Kdprime_0',
                            'R']

        self.param_unit = [ 'J', 'J/K','bar','1/K','K','','K',
                            'J/K-m', 'J/K^2-m', 'J-K/m', 'J-K^(1/2)/m',
                            'J/bar-m','1/K',
                            'bar','', '',
                            'J/K-m']
        
        symdict = dict((p, sym.symbols(p,real=True)) for p in self.param_strs)

        symdict = dict((p, sym.symbols(p,real=True)) for p in self.param_strs)
        
        required_params = self.param_strs
        required_params.remove('R')
        self.param_vals = dict((p, kwargs[p]) for p in required_params)
        self.param_vals['R'] = kwargs.get('R', 8.31446261815324)

        self.syms = types.SimpleNamespace(**symdict)
        self.param_syms = [ symdict[p] for p in self.param_strs ]

    def tait_constants(self):
        """
        returns parameters for the modified Tait equation of state
        derived from K_T and its two first pressure derivatives
        EQ 4 from Holland and Powell, 2011
        """
        Kprime_0 = self.syms.Kprime_0
        K_0 = self.syms.K_0
        Kdprime_0 = self.syms.Kdprime_0
        a = (1.0 + Kprime_0) / (
            1.0 + Kprime_0 + K_0 * Kdprime_0
        )
        b = Kprime_0 / K_0 - Kdprime_0 / (
            1.0 + Kprime_0
        )
        c = (1.0 + Kprime_0 + K_0 * Kdprime_0) / (
            Kprime_0 * Kprime_0
            + Kprime_0
            - K_0 * Kdprime_0
        )
        return a, b, c

    def thermal_energy(self):
        """
        calculate the thermal energy of a substance.  Takes the temperature,
        the Einstein temperature, and n, the number of atoms per molecule.
        Returns thermal energy in J/mol
        """

        T = self.T 
        T_einstein = self.syms.T_einstein
        n = self.syms.n
        R = self.syms.R

        x = T_einstein / T
        E_th = 3.0 * n * R * T_einstein * (1.0 / (sym.exp(x) - 1.0))
        return E_th

    def thermal_energy0(self):
        """
        calculate the thermal energy of a substance.  Takes the temperature,
        the Einstein temperature, and n, the number of atoms per molecule.
        Returns thermal energy in J/mol
        """

        T = self.syms.T_0
        T_einstein = self.syms.T_einstein
        n = self.syms.n
        R = self.syms.R

        x = T_einstein / T
        E_th = 3.0 * n * R * T_einstein * (1.0 / (sym.exp(x) - 1.0))
        return E_th

    def molar_heat_capacity_v(self):
        """
        Heat capacity at constant volume.  In J/K/mol
        """
        T = self.T
        T_einstein = self.syms.T_einstein
        n = self.syms.n
        R = self.syms.R

        x = T_einstein / T
        C_v = (
            3.0
            * n
            * R
            * (x * x * sym.exp(x) / (sym.exp(x) - 1.0)**2.0)
        )
        return C_v

    def thermal_pressure(self):
        """
        Returns thermal pressure [Pa] as a function of T [K]
        EQ 12 - 1 of :cite:`HP2011`.
        """

        # This is basically the mie-gruneisen equation of state for thermal
        # pressure using an Einstein model for heat capacity.  The additional
        # assumption that they make is that alpha*K/Cv, (or gamma / V) is
        # constant over a wide range of compressions.

        # Note that the xi function in HP2011 is just the Einstein heat capacity
        # divided by 3nR. This function is *not* used to calculate the
        # heat capacity - Holland and Powell (2011) prefer the additional
        # freedom provided by their polynomial expression.

        a_0 = self.syms.a_0
        K_0 = self.syms.K_0
        E_th = self.thermal_energy()
        C_V_0 = self.molar_heat_capacity_v()
        P_th = a_0 * K_0 / C_V_0 * E_th
        return P_th

    def thermal_pressure0(self):
        """
        Returns thermal pressure [Pa] as a function of T [K]
        EQ 12 - 1 of :cite:`HP2011`.
        """

        # This is basically the mie-gruneisen equation of state for thermal
        # pressure using an Einstein model for heat capacity.  The additional
        # assumption that they make is that alpha*K/Cv, (or gamma / V) is
        # constant over a wide range of compressions.

        # Note that the xi function in HP2011 is just the Einstein heat capacity
        # divided by 3nR. This function is *not* used to calculate the
        # heat capacity - Holland and Powell (2011) prefer the additional
        # freedom provided by their polynomial expression.

        a_0 = self.syms.a_0
        K_0 = self.syms.K_0
        E_th = self.thermal_energy0()
        C_V_0 = self.molar_heat_capacity_v()
        P_th = a_0 * K_0 / C_V_0 * E_th
        return P_th

    def relative_thermal_pressure(self):
        return self.thermal_pressure() - self.thermal_pressure0()

    def intCpdT(self):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        Cp_0 = self.syms.Cp_0
        Cp_1 = self.syms.Cp_1
        Cp_2 = self.syms.Cp_2
        Cp_3 = self.syms.Cp_3
        T_0 = self.syms.T_0

        T = self.T

        return (
            Cp_0 * T
            + 0.5 * Cp_1 * T**2.0
            - Cp_2 / T
            + 2.0 * Cp_3 * sym.sqrt(T)
        ) - (
            Cp_0 * T_0
            + 0.5 * Cp_1 * T_0 * T_0
            - Cp_2 / T_0
            + 2.0 * Cp_3 * sym.sqrt(T_0)
        )

    def G_PrToP_default(self):
        P = self.P
        Pr = self.syms.Pr
        V_0 = self.syms.V_0
        a,b,c = self.tait_constants()

        theta = self.syms.T_einstein
        n = self.syms.n
        R = self.syms.R
        a_0 = self.syms.a_0
        K_0 = self.syms.K_0
        T = self.T

        CV_0 = 3*n*R*(
            ((theta/T)**2 * sym.exp(theta/T))/(
                (sym.exp(theta/T)-1)**2
            )
        )

        Eth = 3*n*R*(
            0.5 + 1/(sym.exp(theta/T) - 1)
        )

        Pth = a_0*K_0*Eth/CV_0

        V = P*V_0*(1-a+(a*((1-b*Pth)**(1-c)-(1+b*(P-Pth))**(1-c))/(b*(c-1)*P)))
        VdP = sym.integrate(V, (P, Pr, P))
        return VdP

    def intCpoverTdT(self):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        Cp_0 = self.syms.Cp_0
        Cp_1 = self.syms.Cp_1
        Cp_2 = self.syms.Cp_2
        Cp_3 = self.syms.Cp_3
        T = self.T
        T_0 = self.syms.T_0

        return (
            Cp_0 * sym.log(T)
            + Cp_1 * T
            - 0.5 * Cp_2 / T**2.0
            - 2.0 * Cp_3 / sym.sqrt(T)
        ) - (
            Cp_0 * sym.log(T_0)
            + Cp_1 * T_0
            - 0.5 * Cp_2 / (T_0 * T_0)
            - 2.0 * Cp_3 / sym.sqrt(T_0)
        )

    def gibbs_free_energy(self):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # Calculate temperature and pressure integrals
        a, b, c = self.tait_constants()
        P = self.P
        T = self.T
        P_0 = self.syms.P_0
        V_0 = self.syms.V_0
        H_0 = self.syms.H_0
        S_0 = self.syms.S_0
        Pth = self.relative_thermal_pressure()

        psubpth = P - P_0 - Pth

        # EQ 13
        intVdP = (
            (P - P_0)
            * V_0
            * (
                1.0
                - a
                + (
                    a
                    * (
                        (1.0 - b * Pth)**(1.0 - c)
                        - (1.0 + b * (psubpth))**(1.0 - c)
                    )
                    / (b * (c - 1.0) * (P - P_0))
                )
            )
        )
  
        return (
            H_0
            + self.intCpdT()
            - T * (S_0 + self.intCpoverTdT())
            + intVdP
        )

    def add_potential_to_model(self):
        if self.G is None: self.G = self.gibbs_free_energy()
        self.model.add_potential_to_model('G', self.G, self.params())