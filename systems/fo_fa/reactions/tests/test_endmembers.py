import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_fa_binary as db
import pytest
class TestFayalite_berman:
    phase = db.Fayalite_berman()

    def test_name(self):
        test = self.phase.name()
        answer = 'Fayalite_berman'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 1441107.9021681508
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 4.6287570459703536e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 6.93910565957967e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 202.11523159773122
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 176.54853933629937
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -3.3799e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.11889131270454778
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 0.0002254575259
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 7.3154e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 6.18769411314164e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 0.0137005127811399
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 4.870800598538133
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -452.78038873788387
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
        answer = -1988607.1778397025
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
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 452.78038873788387
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 4.870800598538133
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestFayalite_xmelts:
    phase = db.Fayalite_xmelts()

    def test_name(self):
        test = self.phase.name()
        answer = 'Fayalite_xmelts'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0)
        answer = 202174.9855470429
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0)
        answer = 10.03450040340515
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0)
        answer = 0.00010585934253421766
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0)
        answer = 4.946210320204603e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0)
        answer = 240.2
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0)
        answer = 219.36997701235802
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0)
        answer = -2.6750235e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0)
        answer = -0.14129411764705882
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0)
        answer = 0.0005725115
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0)
        answer = 1.46e-09
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0)
        answer = 8.311418685121106e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0)
        answer = -1.15e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0)
        answer = 2.7755575615628914e-17
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0)
        answer = 5.408228374505001
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0)
        answer = -517.7923958463807
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0)
        answer = -2001463.6538405502
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
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0)
        answer = 517.7923958463807
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0)
        answer = 5.408228374505001
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
