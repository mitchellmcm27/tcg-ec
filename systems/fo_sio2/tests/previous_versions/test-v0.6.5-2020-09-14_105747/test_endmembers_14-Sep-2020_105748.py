import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_sio2_poly_linear_rxns as db
import pytest
class TestB_Quartz_berman:
    phase = db.B_Quartz_berman()

    def test_name(self):
        test = self.phase.name()
        answer = 'B_Quartz_berman'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 806559.7498990551
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 1.2398337508475407e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 73.05726934561444
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 73.05726934561444
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -2.934771e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.04297486432094967
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 2.3525676865170618e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.002981213650159617
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 2.367068163771
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -154.12888659874318
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'SiO2'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
        answer = -1074760.335768867
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 60.0843
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 154.12888659874318
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 2.367068163771
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestBeta_Cristobalite_berman:
    phase = db.Beta_Cristobalite_berman()

    def test_name(self):
        test = self.phase.name()
        answer = 'Beta_Cristobalite_berman'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 912405.3648449577
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 3.178481588845701e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 1.0960040772775672e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 73.63336071736431
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 73.59043665134064
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -3.002181e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.04331374159844959
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 8.706516e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 2.3377629288485826e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.003571771808023684
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 2.7392060506356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -156.29266272385198
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'SiO2'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 60.0843
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 156.29266272385198
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 2.7392060506356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestForsterite_berman:
    phase = db.Forsterite_berman()

    def test_name(self):
        test = self.phase.name()
        answer = 'Forsterite_berman'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 1337044.4984269554
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 5.1433326723180724e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 7.479182638846417e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 189.94801785955553
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 162.17297279738264
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -3.4548158e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.11173412815267972
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 0.00023758300657999998
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 7.77148e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 5.734953823236728e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.01423991315765534
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 4.619242458468521
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -373.01618458130076
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'Mg2SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 373.01618458130076
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 4.619242458468521
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestForsterite_xmelts:
    phase = db.Forsterite_xmelts()

    def test_name(self):
        test = self.phase.name()
        answer = 'Forsterite_xmelts'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 375585.4043014601
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = 10.725591783568605
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 0.00010390413350497564
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 2.6625102800782946e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 271.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 236.667434121986
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -1.3260939e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.15941176470588236
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 0.0005175065000000001
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 4.1399999999999997e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 9.377162629757785e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = -6.5e-09
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 2.7755575615628914e-17
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 4.980615135732
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -411.21265315352593
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'Mg2SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
        answer = -2543709.4237597957
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 411.21265315352593
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 4.980615135732
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestOrthoenstatite_berman:
    phase = db.Orthoenstatite_berman()

    def test_name(self):
        test = self.phase.name()
        answer = 'Orthoenstatite_berman'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 1399378.163884626
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 4.3568139995068164e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 7.146031185909278e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 136.65231010553666
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 121.8177977885435
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -2.3475568999999997e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.08038371182678627
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 0.00014312656214999998
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 4.6995e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 4.17896152261526e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.00934136594232686
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 3.285119864336684
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -266.2955829601318
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'MgSiO3'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 100.3887
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 266.2955829601318
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 3.285119864336684
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestQuartz4_liquid:
    phase = db.Quartz4_liquid()

    def test_name(self):
        test = self.phase.name()
        answer = 'Quartz4_liquid'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 291926.44334738934
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = 4.722229324452722
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 1.2111410674701189e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 3.425520444579975e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 162.746
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 162.74209702353176
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -1.8365835e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.09573294117647059
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 6.493500000000001e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 3.6e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 5.6313494809688585e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 6.5e-09
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 5.3614728906549995
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -314.4903025207027
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'Si2O4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
        answer = -2150052.4704855355
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 120.1686
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 314.4903025207027
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 5.3614728906549995
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestQuartz4_xmelts:
    phase = db.Quartz4_xmelts()

    def test_name(self):
        test = self.phase.name()
        answer = 'Quartz4_xmelts'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 291926.44334738934
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = 4.722229324452722
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 1.2111410674701189e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 3.425520444579975e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 165.2
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 165.19609702353173
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -1.8365835e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.09717647058823528
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 6.493500000000001e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 3.6e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 5.716262975778546e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 6.5e-09
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 5.3614728906549995
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -318.65618011721415
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'Si2O4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
        answer = -2148722.9341392373
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 120.1686
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 318.65618011721415
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 5.3614728906549995
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestSilica_polymorph_berman:
    phase = db.Silica_polymorph_berman()

    def test_name(self):
        test = self.phase.name()
        answer = 'Silica_polymorph_berman'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 912405.3648449577
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 3.178481588845701e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 1.0960040772775672e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 73.63336071736431
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 73.59043665134064
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -3.002181e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.04331374159844959
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 8.706516e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 2.3377629288485826e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.003571771808023684
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 2.7392060506356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -156.29266272385198
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'SiO2'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
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

    def test_mw(self):
        test = self.phase.mw()
        answer = 60.0843
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 156.29266272385198
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 2.7392060506356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

