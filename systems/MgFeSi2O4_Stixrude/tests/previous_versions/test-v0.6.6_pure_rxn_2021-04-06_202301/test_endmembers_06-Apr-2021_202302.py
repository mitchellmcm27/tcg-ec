import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_Mg2SiO4_stixrude as db
import pytest
class TestFayalite_stixrude:
    phase = db.Fayalite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'Fayalite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 1427805.4194400315
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.692110212623921
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 2.5361303641362633e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 7.003755458444669e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 180.34855139356458
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 173.3703477497622
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -3.1304986318008575e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.10608738317268504
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 0.00011335850719094473
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 1.2480092168068203e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 3.2530087321578803e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 5.7618459376512696e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -6.620961655580243e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.00813600223261346
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 4.469742912034868
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -403.8210368856753
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
        test = self.phase.g(1700.0, 100000.0)
        answer = -1335279.3432542302
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
        test = self.phase.s(1700.0, 100000.0)
        answer = 403.8210368856753
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 4.469742912034868
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestFePerovskite_stixrude:
    phase = db.FePerovskite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'FePerovskite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 2689740.4683418437
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.274109553850979
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 2.8105277314066333e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 3.7178308159094415e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 132.30144926420007
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 123.09622165696729
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -9.475215712864337e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.07782438192011769
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 7.162874762376963e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 1.8579237032048664e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 1.2322268205765097e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 4.036198698626725e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -1.4445585010429686e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.009209004043463365
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 2.548587114915972
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -250.54072038317454
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'FeSiO3'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0)
        answer = -1032107.9475255425
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
        answer = 131.9307
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0)
        answer = 250.54072038317454
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 2.548587114915972
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestFeRingwoodite_stixrude:
    phase = db.FeRingwoodite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'FeRingwoodite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 2161744.2017083573
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.585612787379432
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 2.3891603140594098e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 4.6258942163912456e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 181.8370444223005
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 173.17617410756017
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -1.909912294822873e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.10696296730723558
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 9.864226125094362e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 4.93491807601715e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 2.0425604604346202e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 5.762412880581958e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -3.159702976233127e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.009001948337342289
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 4.1287418291048485
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -391.08671430995014
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
        test = self.phase.g(1700.0, 100000.0)
        answer = -1347469.3172791
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
        test = self.phase.s(1700.0, 100000.0)
        answer = 391.08671430995014
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 4.1287418291048485
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestFeWadsleyite_stixrude:
    phase = db.FeWadsleyite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'FeWadsleyite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 1770679.3920615853
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.590644530028833
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 2.7032540278350814e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 5.647549773738031e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 182.43699126876467
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 173.209252249792
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -2.369151310574721e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.1073158772169204
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 0.0001134017065709345
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 7.480226445654897e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 2.402442300559525e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 5.7582905104922745e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -4.188221815674019e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.009424938538551739
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 4.195007402310355
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -393.12840946831494
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
        test = self.phase.g(1700.0, 100000.0)
        answer = -1344878.4982517625
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
        test = self.phase.s(1700.0, 100000.0)
        answer = 393.12840946831494
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 4.195007402310355
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestForsterite_stixrude:
    phase = db.Forsterite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'Forsterite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 1429530.0272116363
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.371247101982014
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 2.6450925480386394e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 6.995306016415379e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 179.65452924987395
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 172.50308587748762
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -2.942229018403993e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.10567913485286703
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 0.00011125271765010549
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 1.1054989253562702e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 2.344995657194957e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 5.71274705310251e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -4.981947378765701e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.008562434950124362
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 4.206004728741926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -357.46364178007343
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
        test = self.phase.g(1700.0, 100000.0)
        answer = -1984446.1221989095
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
        test = self.phase.s(1700.0, 100000.0)
        answer = 357.46364178007343
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 4.206004728741926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestMgPerovskite_stixrude:
    phase = db.MgPerovskite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'MgPerovskite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 2455273.6800218695
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.302105774495801
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 3.207220234066898e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 4.0728657181349043e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 133.48390635877956
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 122.9723459619135
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -9.971509001699658e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.07851994491692915
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 7.852167895453328e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 2.1533239202025615e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 1.544309042514691e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 3.990400107339729e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -1.7460043586313647e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.010683143092153763
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 2.4482783601974316
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -245.96599745551111
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
        test = self.phase.g(1700.0, 100000.0)
        answer = -1363449.0604146551
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
        test = self.phase.s(1700.0, 100000.0)
        answer = 245.96599745551111
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 2.4482783601974316
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestMgRingwoodite_stixrude:
    phase = db.MgRingwoodite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'MgRingwoodite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 1916858.2072982239
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.5245310300166635
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 2.4539534710325686e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 5.216869960399844e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 179.7983142829723
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 172.1994273285761
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -2.0201724176636655e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.10576371428410136
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 9.50265035172556e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 5.822290435919738e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 2.0705894836176042e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 5.666994947427094e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -3.4060101483156335e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.009424800177840767
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 3.8723840789560926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -345.72095455231477
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
        test = self.phase.g(1700.0, 100000.0)
        answer = -1966692.8925920525
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
        test = self.phase.s(1700.0, 100000.0)
        answer = 345.72095455231477
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 3.8723840789560926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestMgWadsleyite_stixrude:
    phase = db.MgWadsleyite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'MgWadsleyite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 1757071.3405069287
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.6238090125529325
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 2.8670619149850206e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 5.6912885489982e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 182.13818429275685
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 172.37672081299118
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -2.262627591551502e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.10714010840750403
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 0.00011398285888480343
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 7.2419287413493e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 2.64020860816539e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 5.645682346622832e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -4.2741988691355865e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.011163508514915876
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 3.975598095355361
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -352.3823880256589
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
        test = self.phase.g(1700.0, 100000.0)
        answer = -1975655.3236001502
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
        test = self.phase.s(1700.0, 100000.0)
        answer = 352.3823880256589
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 3.975598095355361
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestPericlase_stixrude:
    phase = db.Periclase_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'Periclase_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 1629319.3004791553
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 4.1904091850190825
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 3.636368711886252e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 6.137532402064573e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 53.43230243192754
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 49.365418571070265
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -6.814980540518679e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.031430766136427965
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 4.0377435728599276e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 2.1710009562172374e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 1.0148203879193656e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 1.6002932955665014e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -1.4575289337765742e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.004225780111797441
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 1.110377932705695
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -105.63270649432897
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'MgO'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0)
        answer = -563891.0645437283
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
        answer = 40.3044
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0)
        answer = 105.63270649432897
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 1.110377932705695
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))


class TestWuestite_stixrude:
    phase = db.Wuestite_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'Wuestite_stixrude'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0)
        answer = 1805805.4394082627
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0)
        answer = 5.32196804073994
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0)
        answer = 3.415191676976622e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0)
        answer = 5.537695136900718e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0)
        answer = 54.05994713620336
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0)
        answer = 49.704773175038255
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0)
        answer = -6.735727437930723e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0)
        answer = -0.031799968903649035
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0)
        answer = 4.1540387680638215e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0)
        answer = 2.358119688003987e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0)
        answer = 1.110575181748639e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0)
        answer = 1.6305900476593658e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0)
        answer = -1.613064198088169e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0)
        answer = 0.004079938093439817
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0)
        answer = 1.216341324578678
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0)
        answer = -131.82994540737673
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_elements(self):
        test = self.phase.elements()
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_formula(self):
        test = self.phase.formula()
        answer = 'FeO'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0)
        answer = -261701.67654927325
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
        answer = 71.8464
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0)
        answer = 131.82994540737673
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0)
        answer = 1.216341324578678
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

