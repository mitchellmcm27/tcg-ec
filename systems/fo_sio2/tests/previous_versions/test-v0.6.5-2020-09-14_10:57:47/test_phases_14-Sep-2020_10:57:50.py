import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_sio2_poly_linear_rxns as db
import pytest
class TestLiquid:
    phase = db.Liquid()

    def test_name(self):
        test = self.phase.name()
        answer = 'Liquid'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 330737.88899648655
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 7.977299588566115
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([0.28071311, 0.71928689]),)
        answer = 134.93160385355992
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Lq'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 7.696919774517023e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 3.0235423072759076e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([0.25, 0.75]),)
        answer =[0.2807131061141594, 0.7192868938858406]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 4.        , 0.        ,
       0.        , 0.        , 1.43857379, 0.        , 1.28071311,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        ]),)
        answer =[0.28071310611415934, 0.7192868938858406]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 4.        , 0.        ,
       0.        , 0.        , 1.43857379, 0.        , 1.28071311,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        ]),)
        answer = 134.93160385355992
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 4.        , 0.        ,
       0.        , 0.        , 1.43857379, 0.        , 1.28071311,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        , 0.        , 0.        , 0.        , 0.        ,
       0.        ]),)
        answer = 0.9999999999999999
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([0.28071311, 0.71928689]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.438573787771681, 0.0, 1.2807131061141592, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_mole_frac(self):
        test = self.phase.conv_moles_to_mole_frac(array([0.28071311, 0.71928689]),)
        answer =[0.2807131061141594, 0.7192868938858406]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_tot_moles(self):
        test = self.phase.conv_moles_to_tot_moles(array([0.28071311, 0.71928689]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 240.61168341071777
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 224.42382841970604
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[112037.45259824768, -43724.390903418884, 17064.13628245145]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[5.166457164616439, 4.740186238641632]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[-325.0531796293239, -413.9522268719305]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = -1.4693950212549747e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = -0.14153628435924576
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 0.0003740584535052851
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[2.380385431413572, -0.9289831274293263, 0.36255038350486757]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[21.304612647865667, -8.31446261815324, 3.2448507640711894]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[-603351.5583045483, 79705.66112014651, 29682.130144409835, -35307.60736398399]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[-1.8365835e-05, -1.3260939e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[-0.09573294117647059, -0.15941176470588236]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[6.493500000000001e-06, 0.0005175065000000001]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 3.988414922698353e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 8.325663785837986e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = -2.8507296205159275e-09
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[-2170810.1603127318, -2556227.53702241]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 4.859846074318178
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = -388.997099189869
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[5.166457164616439, 4.740186238641632]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[-325.0531796293239, -413.9522268719305]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 1000.0, array([0.25, 0.75]))
        answer =[[125801.52527988631, -41933.841759962095], [-49096.03833379451, 16365.346111264835]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 1000.0, array([0.25, 0.75]))
        answer = 8.394739392997513e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 1000.0, array([0.25, 0.75]))
        answer = -0.002137017745059897
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 1000.0, array([0.25, 0.75]))
        answer =[-5.377832467156563, 1.7926108223855233]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 1000.0, array([0.25, 0.75]))
        answer =[-71.79958510781829, 23.933195035939445]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 1000.0, array([0.25, 0.75]))
        answer =[0.344279006106694, -0.11475966870223148]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'Si2O4'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 120.1686
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'Quartz4_liquid'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 2
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 'Mg1.439Si1.281O4.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = -2448035.828055864
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_names(self):
        test = self.phase.get_param_names()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_number(self):
        test = self.phase.get_param_number()
        answer = 0
        assert(test == answer)

    def test_get_param_units(self):
        test = self.phase.get_param_units()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_values(self):
        test = self.phase.get_param_values()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_mu(self):
        test = self.phase.mu(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer =[-2170810.1603127318, -2556227.53702241]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, array([0.25, 0.75]))
        answer = 27.764583855156445
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 388.997099189869
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'Si2O4'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 120.1686
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'Quartz4_liquid'
        assert(test == answer)

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 2
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([0.28071311, 0.71928689]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0, array([0.28071311, 0.71928689]))
        answer = 4.859846074318178
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([0.28071311, 0.71928689]),)
        answer =[0.24999999999999994, 0.75]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestOlivine:
    phase = db.Olivine()

    def test_name(self):
        test = self.phase.name()
        answer = 'Olivine'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0, array([1.]))
        answer = 1337044.4984269554
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0, array([1.]))
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Ol'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, array([1.]))
        answer = 5.1433326723180724e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0, array([1.]))
        answer = 7.479182638846417e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 2., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 2., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 2., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_mole_frac(self):
        test = self.phase.conv_moles_to_mole_frac(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_tot_moles(self):
        test = self.phase.conv_moles_to_tot_moles(array([1.]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0, array([1.]))
        answer = 189.94801785955553
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0, array([1.]))
        answer = 162.17297279738264
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 1000.0, array([1.]))
        answer =[4.619242458468521]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 1000.0, array([1.]))
        answer =[-373.01618458130076]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0, array([1.]))
        answer = -3.4548158e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0, array([1.]))
        answer = -0.11173412815267972
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0, array([1.]))
        answer = 0.00023758300657999998
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 1000.0, array([1.]))
        answer =[-3.4548158e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 1000.0, array([1.]))
        answer =[-0.11173412815267972]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 1000.0, array([1.]))
        answer =[0.00023758300657999998]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0, array([1.]))
        answer = 7.77148e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0, array([1.]))
        answer = 5.734953823236728e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0, array([1.]))
        answer = 0.01423991315765534
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 1000.0, array([1.]))
        answer =[-2566345.4730803007]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0, array([1.]))
        answer = 4.619242458468521
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0, array([1.]))
        answer = -373.01618458130076
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 1000.0, array([1.]))
        answer =[4.619242458468521]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 1000.0, array([1.]))
        answer =[-373.01618458130076]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 1000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 1000.0, array([1.]))
        answer = 2.2780129001376464e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 1000.0, array([1.]))
        answer = -0.0015665586392268071
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'Mg2SiO4'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'Forsterite_berman'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 1000.0, array([1.]))
        answer = 'Mg2.000Si1.000O4.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0, array([1.]))
        answer = -2566345.4730803007
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_names(self):
        test = self.phase.get_param_names()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_number(self):
        test = self.phase.get_param_number()
        answer = 0
        assert(test == answer)

    def test_get_param_units(self):
        test = self.phase.get_param_units()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_values(self):
        test = self.phase.get_param_values()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_mu(self):
        test = self.phase.mu(1700.0, 1000.0, array([1.]))
        answer =[-2566345.4730803007]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, array([1.]))
        answer = 30.458046154746736
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, array([1.]))
        answer = 373.01618458130076
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'Mg2SiO4'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'Forsterite_berman'
        assert(test == answer)

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 1
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([1.]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0, array([1.]))
        answer = 4.619242458468521
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestOrthopyroxene:
    phase = db.Orthopyroxene()

    def test_name(self):
        test = self.phase.name()
        answer = 'Orthopyroxene'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0, array([1.]))
        answer = 1399378.163884626
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0, array([1.]))
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 100.3887
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Opx'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, array([1.]))
        answer = 4.3568139995068164e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0, array([1.]))
        answer = 7.146031185909278e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0., 0., 0., 0., 0., 0., 0., 0., 3., 0., 0., 0., 1., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0., 0., 0., 0., 0., 0., 0., 0., 3., 0., 0., 0., 1., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 100.3887
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0., 0., 0., 0., 0., 0., 0., 0., 3., 0., 0., 0., 1., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_mole_frac(self):
        test = self.phase.conv_moles_to_mole_frac(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_tot_moles(self):
        test = self.phase.conv_moles_to_tot_moles(array([1.]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0, array([1.]))
        answer = 136.65231010553666
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0, array([1.]))
        answer = 121.8177977885435
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 1000.0, array([1.]))
        answer =[3.285119864336684]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 1000.0, array([1.]))
        answer =[-266.2955829601318]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0, array([1.]))
        answer = -2.3475568999999997e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0, array([1.]))
        answer = -0.08038371182678627
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0, array([1.]))
        answer = 0.00014312656214999998
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 1000.0, array([1.]))
        answer =[-2.3475568999999997e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 1000.0, array([1.]))
        answer =[-0.08038371182678627]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 1000.0, array([1.]))
        answer =[0.00014312656214999998]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0, array([1.]))
        answer = 4.6995e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0, array([1.]))
        answer = 4.17896152261526e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0, array([1.]))
        answer = 0.00934136594232686
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 1000.0, array([1.]))
        answer =[-1823995.2293042827]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0, array([1.]))
        answer = 3.285119864336684
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0, array([1.]))
        answer = -266.2955829601318
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 1000.0, array([1.]))
        answer =[3.285119864336684]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 1000.0, array([1.]))
        answer =[-266.2955829601318]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 1000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 1000.0, array([1.]))
        answer = 2.1837278715482756e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 1000.0, array([1.]))
        answer = -0.0013313818418026052
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'MgSiO3'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 100.3887
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'Orthoenstatite_berman'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 1000.0, array([1.]))
        answer = 'Mg1.000Si1.000O3.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0, array([1.]))
        answer = -1823995.2293042827
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_names(self):
        test = self.phase.get_param_names()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_number(self):
        test = self.phase.get_param_number()
        answer = 0
        assert(test == answer)

    def test_get_param_units(self):
        test = self.phase.get_param_units()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_values(self):
        test = self.phase.get_param_values()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_mu(self):
        test = self.phase.mu(1700.0, 1000.0, array([1.]))
        answer =[-1823995.2293042827]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, array([1.]))
        answer = 30.55861099310908
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, array([1.]))
        answer = 266.2955829601318
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'MgSiO3'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 100.3887
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'Orthoenstatite_berman'
        assert(test == answer)

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 1
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([1.]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0, array([1.]))
        answer = 3.285119864336684
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestSilica_polymorph:
    phase = db.Silica_polymorph()

    def test_name(self):
        test = self.phase.name()
        answer = 'Silica_polymorph'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0, array([1.]))
        answer = 912405.3648449577
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0, array([1.]))
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 60.0843
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'pQz'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, array([1.]))
        answer = 3.178481588845701e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0, array([1.]))
        answer = 1.0960040772775672e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 60.0843
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_mole_frac(self):
        test = self.phase.conv_moles_to_mole_frac(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_tot_moles(self):
        test = self.phase.conv_moles_to_tot_moles(array([1.]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0, array([1.]))
        answer = 73.63336071736431
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0, array([1.]))
        answer = 73.59043665134064
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 1000.0, array([1.]))
        answer =[2.7392060506356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 1000.0, array([1.]))
        answer =[-156.29266272385198]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0, array([1.]))
        answer = -3.002181e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0, array([1.]))
        answer = -0.04331374159844959
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0, array([1.]))
        answer = 8.706516e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 1000.0, array([1.]))
        answer =[-3.002181e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 1000.0, array([1.]))
        answer =[-0.04331374159844959]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 1000.0, array([1.]))
        answer =[8.706516e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0, array([1.]))
        answer = 2.3377629288485826e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0, array([1.]))
        answer = 0.003571771808023684
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 1000.0, array([1.]))
        answer =[-1075500.4686895711]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0, array([1.]))
        answer = 2.7392060506356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0, array([1.]))
        answer = -156.29266272385198
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 1000.0, array([1.]))
        answer =[2.7392060506356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 1000.0, array([1.]))
        answer =[-156.29266272385198]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 1000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 1000.0, array([1.]))
        answer = 2.404077552511547e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 1000.0, array([1.]))
        answer = -6.97197793077187e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 1000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'SiO2'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 60.0843
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'Silica_polymorph_berman'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 1000.0, array([1.]))
        answer = 'Si1.000O2.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0, array([1.]))
        answer = -1075500.4686895711
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_names(self):
        test = self.phase.get_param_names()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_number(self):
        test = self.phase.get_param_number()
        answer = 0
        assert(test == answer)

    def test_get_param_units(self):
        test = self.phase.get_param_units()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_param_values(self):
        test = self.phase.get_param_values()
        answer =[]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_mu(self):
        test = self.phase.mu(1700.0, 1000.0, array([1.]))
        answer =[-1075500.4686895711]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, array([1.]))
        answer = 21.93493256414871
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, array([1.]))
        answer = 156.29266272385198
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'SiO2'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 60.0843
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'Silica_polymorph_berman'
        assert(test == answer)

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 1
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([1.]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0, array([1.]))
        answer = 2.7392060506356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

