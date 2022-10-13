import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_h2o_hydration as db
import pytest
class TestBrucite:
    phase = db.Brucite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Brucite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(423.15, 10.0, array([1.]))
        answer = 496393.9624777537
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(423.15, 10.0, array([1.]))
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 58.319599999999994
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Brc'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(423.15, 10.0, array([1.]))
        answer = 3.299333900718681e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(423.15, 10.0, array([1.]))
        answer = 2.01452893385023e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0., 2., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 1., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0., 2., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 1., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 58.319599999999994
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0., 2., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 1., 0., 0., 0., 0.,
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
        answer =[0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(423.15, 10.0, array([1.]))
        answer = 93.66123573296828
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(423.15, 10.0, array([1.]))
        answer = 93.09460781874695
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(423.15, 10.0, array([1.]))
        answer =[2.4781328856164]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(423.15, 10.0, array([1.]))
        answer =[-93.14855222831395]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(423.15, 10.0, array([1.]))
        answer = -4.9922704e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(423.15, 10.0, array([1.]))
        answer = -0.22134287069116929
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(423.15, 10.0, array([1.]))
        answer = 8.176187839999999e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(423.15, 10.0, array([1.]))
        answer =[-4.9922704e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(423.15, 10.0, array([1.]))
        answer =[-0.22134287069116929]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(423.15, 10.0, array([1.]))
        answer =[8.176187839999999e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(423.15, 10.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(423.15, 10.0, array([1.]))
        answer = 5.4295999999999995e-09
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(423.15, 10.0, array([1.]))
        answer = 0.00030028432933008306
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(423.15, 10.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(423.15, 10.0, array([1.]))
        answer = 0.09427755673514465
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(423.15, 10.0, array([1.]))
        answer =[-954530.8351882802]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(423.15, 10.0, array([1.]))
        answer = 2.4781328856164
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(423.15, 10.0, array([1.]))
        answer = -93.14855222831395
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(423.15, 10.0, array([1.]))
        answer =[2.4781328856164]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(423.15, 10.0, array([1.]))
        answer =[-93.14855222831395]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(423.15, 10.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(423.15, 10.0, array([1.]))
        answer = 4.740929039459025e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(423.15, 10.0, array([1.]))
        answer = -0.0007764548643584644
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(423.15, 10.0, array([1.]))
        answer = 'Mg1.000O2.000H2.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(423.15, 10.0, array([1.]))
        answer = -954530.8351882802
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
        test = self.phase.mu(423.15, 10.0, array([1.]))
        answer =[-954530.8351882802]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(423.15, 10.0, array([1.]))
        answer = 23.533685517229163
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(423.15, 10.0, array([1.]))
        answer = 93.14855222831395
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 1
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([1.]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(423.15, 10.0, array([1.]))
        answer = 2.4781328856164
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestChrysotile:
    phase = db.Chrysotile()

    def test_name(self):
        test = self.phase.name()
        answer = 'Chrysotile'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(423.15, 10.0, array([1.]))
        answer = 554318.238579241
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(423.15, 10.0, array([1.]))
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 277.11220000000003
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Ctl'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(423.15, 10.0, array([1.]))
        answer = 2.8726266313116093e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(423.15, 10.0, array([1.]))
        answer = 1.8040178554526272e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0., 4., 0., 0., 0., 0., 0., 0., 9., 0., 0., 0., 3., 0., 2., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0., 4., 0., 0., 0., 0., 0., 0., 9., 0., 0., 0., 3., 0., 2., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 277.11220000000003
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0., 4., 0., 0., 0., 0., 0., 0., 9., 0., 0., 0., 3., 0., 2., 0., 0.,
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
        answer =[0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 3.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(423.15, 10.0, array([1.]))
        answer = 330.911783074336
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(423.15, 10.0, array([1.]))
        answer = 328.8296121227136
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(423.15, 10.0, array([1.]))
        answer =[10.757330334256]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(423.15, 10.0, array([1.]))
        answer =[-326.23863576394433]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(423.15, 10.0, array([1.]))
        answer = -1.9406416e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(423.15, 10.0, array([1.]))
        answer = -0.7820200474402362
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(423.15, 10.0, array([1.]))
        answer = 0.000309017936
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(423.15, 10.0, array([1.]))
        answer =[-1.9406416e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(423.15, 10.0, array([1.]))
        answer =[-0.7820200474402362]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(423.15, 10.0, array([1.]))
        answer =[0.000309017936]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(423.15, 10.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(423.15, 10.0, array([1.]))
        answer = 1.43648e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(423.15, 10.0, array([1.]))
        answer = 0.0010178280172080047
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(423.15, 10.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(423.15, 10.0, array([1.]))
        answer = 0.3513261219586691
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(423.15, 10.0, array([1.]))
        answer =[-4463217.900778551]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(423.15, 10.0, array([1.]))
        answer = 10.757330334256
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(423.15, 10.0, array([1.]))
        answer = -326.23863576394433
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(423.15, 10.0, array([1.]))
        answer =[10.757330334256]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(423.15, 10.0, array([1.]))
        answer =[-326.23863576394433]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(423.15, 10.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(423.15, 10.0, array([1.]))
        answer = 4.647206520857805e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(423.15, 10.0, array([1.]))
        answer = -0.000739997621014215
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(423.15, 10.0, array([1.]))
        answer = 'Mg3.000Si2.000O9.000H4.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(423.15, 10.0, array([1.]))
        answer = -4463217.900778551
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
        test = self.phase.mu(423.15, 10.0, array([1.]))
        answer =[-4463217.900778551]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(423.15, 10.0, array([1.]))
        answer = 25.76031332955861
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(423.15, 10.0, array([1.]))
        answer = 326.23863576394433
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 1
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([1.]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(423.15, 10.0, array([1.]))
        answer = 10.757330334256
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestOlivine:
    phase = db.Olivine()

    def test_name(self):
        test = self.phase.name()
        answer = 'Olivine'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(423.15, 10.0, array([1.]))
        answer = 1268564.2655124478
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(423.15, 10.0, array([1.]))
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
        test = self.phase.alpha(423.15, 10.0, array([1.]))
        answer = 3.156827010221178e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(423.15, 10.0, array([1.]))
        answer = 7.88292739426994e-07
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
        test = self.phase.cp(423.15, 10.0, array([1.]))
        answer = 139.8194724211109
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(423.15, 10.0, array([1.]))
        answer = 137.47499660681822
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(423.15, 10.0, array([1.]))
        answer =[4.3826558678078]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(423.15, 10.0, array([1.]))
        answer =[-139.38692215456356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(423.15, 10.0, array([1.]))
        answer = -3.4548158e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(423.15, 10.0, array([1.]))
        answer = -0.3304253158953348
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(423.15, 10.0, array([1.]))
        answer = 0.0001383528642
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(423.15, 10.0, array([1.]))
        answer =[-3.4548158e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(423.15, 10.0, array([1.]))
        answer =[-0.3304253158953348]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(423.15, 10.0, array([1.]))
        answer =[0.0001383528642]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(423.15, 10.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(423.15, 10.0, array([1.]))
        answer = 7.77148e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(423.15, 10.0, array([1.]))
        answer = 0.0004834994635937993
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(423.15, 10.0, array([1.]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(423.15, 10.0, array([1.]))
        answer = 0.12583251787561864
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(423.15, 10.0, array([1.]))
        answer =[-2217084.332036725]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(423.15, 10.0, array([1.]))
        answer = 4.3826558678078
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(423.15, 10.0, array([1.]))
        answer = -139.38692215456356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(423.15, 10.0, array([1.]))
        answer =[4.3826558678078]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(423.15, 10.0, array([1.]))
        answer =[-139.38692215456356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(423.15, 10.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(423.15, 10.0, array([1.]))
        answer = 2.5305968016364413e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(423.15, 10.0, array([1.]))
        answer = -0.0010134123956529342
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(423.15, 10.0, array([1.]))
        answer = 'Mg2.000Si1.000O4.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(423.15, 10.0, array([1.]))
        answer = -2217084.332036725
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
        test = self.phase.mu(423.15, 10.0, array([1.]))
        answer =[-2217084.332036725]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(423.15, 10.0, array([1.]))
        answer = 32.102246729760815
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(423.15, 10.0, array([1.]))
        answer = 139.38692215456356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 1
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([1.]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(423.15, 10.0, array([1.]))
        answer = 4.3826558678078
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestWater:
    phase = db.Water()

    def test_name(self):
        test = self.phase.name()
        answer = 'Water'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(423.15, 10.0, array([1.]))
        answer = 16164.199305489332
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(423.15, 10.0, array([1.]))
        answer = 7.666990737958184
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 18.01528
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'H2O'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(423.15, 10.0, array([1.]))
        answer = 0.0010244809642597453
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(423.15, 10.0, array([1.]))
        answer = 6.186511197374323e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0., 2., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0., 2., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0.]),)
        answer = 18.01528
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0., 2., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
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
        answer =[0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(423.15, 10.0, array([1.]))
        answer = 77.62011079817877
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(423.15, 10.0, array([1.]))
        answer = 63.52125218290239
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(423.15, 10.0, array([1.]))
        answer =[1.9639373384173742]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(423.15, 10.0, array([1.]))
        answer =[-96.50695358089571]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(423.15, 10.0, array([1.]))
        answer = -0.0001214992033506061
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(423.15, 10.0, array([1.]))
        answer = -0.18343403237192196
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(423.15, 10.0, array([1.]))
        answer = 0.0020120164182075495
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(423.15, 10.0, array([1.]))
        answer =[-0.0001214992033506061]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(423.15, 10.0, array([1.]))
        answer =[-0.18343403237192196]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(423.15, 10.0, array([1.]))
        answer =[0.0020120164182075495]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(423.15, 10.0, array([1.]))
        answer = 6.514597167528076e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(423.15, 10.0, array([1.]))
        answer = 1.3777537519832142e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(423.15, 10.0, array([1.]))
        answer = 0.0003222764336543916
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(423.15, 10.0, array([1.]))
        answer = -8.429019059527245e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(423.15, 10.0, array([1.]))
        answer = 0.04706275947106617
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(423.15, 10.0, array([1.]))
        answer =[-247654.0271157971]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(423.15, 10.0, array([1.]))
        answer = 1.9639373384173742
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(423.15, 10.0, array([1.]))
        answer = -96.50695358089571
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(423.15, 10.0, array([1.]))
        answer =[1.9639373384173742]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(423.15, 10.0, array([1.]))
        answer =[-96.50695358089571]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(423.15, 10.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(423.15, 10.0, array([1.]))
        answer = 0.0005674912802138908
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(423.15, 10.0, array([1.]))
        answer = -0.009397607074715634
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(423.15, 10.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(423.15, 10.0, array([1.]))
        answer = 'H2.000O1.000'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(423.15, 10.0, array([1.]))
        answer = -247654.0271157971
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
        test = self.phase.mu(423.15, 10.0, array([1.]))
        answer =[-247654.0271157971]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(423.15, 10.0, array([1.]))
        answer = 9.173042157504625
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(423.15, 10.0, array([1.]))
        answer = 96.50695358089571
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 1
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([1.]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(423.15, 10.0, array([1.]))
        answer = 1.9639373384173742
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

