import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_Mg2SiO4_stixrude as db
import pytest
class TestFayalite:
    phase = db.Fayalite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Fayalite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1427805.4194400315
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.692110212623921
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Fa'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 2.5361303641362633e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 7.003755458444669e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 180.34855139356458
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 173.3703477497622
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[4.469742912034868]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-403.8210368856753]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -3.1304986318008575e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10608738317268504
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.00011335850719094473
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-3.1304986318008575e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10608738317268504]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.00011335850719094473]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 1.2480092168068203e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 3.2530087321578803e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.7618459376512696e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -6.620961655580243e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.00813600223261346
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1335279.3432542302]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.469742912034868
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -403.8210368856753
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.469742912034868]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-403.8210368856753]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 3.193035940810485e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0011562304611169552
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'Fayalite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg0.000SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1335279.3432542302]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 45.590340207560104
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 403.8210368856753
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'Fayalite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 4.469742912034868
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestFePerovskite:
    phase = db.FePerovskite()

    def test_name(self):
        test = self.phase.name()
        answer = 'FePerovskite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 2689740.4683418437
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.274109553850979
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 131.9307
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'FePv'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 2.8105277314066333e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 3.7178308159094415e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 131.9307
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 132.30144926420007
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 123.09622165696729
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[2.548587114915972]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-250.54072038317454]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -9.475215712864337e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.07782438192011769
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 7.162874762376963e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-9.475215712864337e-07]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.07782438192011769]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[7.162874762376963e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 1.8579237032048664e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 1.2322268205765097e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 4.036198698626725e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -1.4445585010429686e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.009209004043463365
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1032107.9475255425]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 2.548587114915972
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -250.54072038317454
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[2.548587114915972]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-250.54072038317454]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 1.924580168964229e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0014549037339306897
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'FeSiO3'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 131.9307
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'FePerovskite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Fe1.000SiO3'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1032107.9475255425]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 51.7662116503127
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 250.54072038317454
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'FeSiO3'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 131.9307
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'FePerovskite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 2.548587114915972
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestFeRingwoodite:
    phase = db.FeRingwoodite()

    def test_name(self):
        test = self.phase.name()
        answer = 'FeRingwoodite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 2161744.2017083573
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.585612787379432
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'FeRi'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 2.3891603140594098e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 4.6258942163912456e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 181.8370444223005
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 173.17617410756017
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[4.1287418291048485]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-391.08671430995014]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -1.909912294822873e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10696296730723558
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 9.864226125094362e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-1.909912294822873e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10696296730723558]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[9.864226125094362e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 4.93491807601715e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.0425604604346202e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.762412880581958e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -3.159702976233127e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.009001948337342289
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1347469.3172791]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.1287418291048485
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -391.08671430995014
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.1287418291048485]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-391.08671430995014]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 2.2831442297455458e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0011791877050827152
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'FeRingwoodite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg0.000SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1347469.3172791]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 49.355738003163275
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 391.08671430995014
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'FeRingwoodite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 4.1287418291048485
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestFeWadsleyite:
    phase = db.FeWadsleyite()

    def test_name(self):
        test = self.phase.name()
        answer = 'FeWadsleyite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1770679.3920615853
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.590644530028833
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'FeWa'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 2.7032540278350814e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 5.647549773738031e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 182.43699126876467
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 173.209252249792
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[4.195007402310355]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-393.12840946831494]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -2.369151310574721e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.1073158772169204
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.0001134017065709345
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-2.369151310574721e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.1073158772169204]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.0001134017065709345]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 7.480226445654897e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.402442300559525e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.7582905104922745e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -4.188221815674019e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.009424938538551739
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1344878.4982517625]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.195007402310355
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -393.12840946831494
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.195007402310355]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-393.12840946831494]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 2.743359438088664e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0013131353857735065
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'FeWadsleyite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg0.000SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1344878.4982517625]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 48.57610022041248
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 393.12840946831494
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'Fe2SiO4'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 203.77710000000002
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'FeWadsleyite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 4.195007402310355
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestForsterite:
    phase = db.Forsterite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Forsterite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1429530.0272116363
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.371247101982014
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Fo'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 2.6450925480386394e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 6.995306016415379e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 179.65452924987395
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 172.50308587748762
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[4.206004728741926]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-357.46364178007343]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -2.942229018403993e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10567913485286703
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.00011125271765010549
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-2.942229018403993e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10567913485286703]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.00011125271765010549]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 1.1054989253562702e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.344995657194957e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.71274705310251e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -4.981947378765701e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.008562434950124362
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1984446.1221989095]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.206004728741926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -357.46364178007343
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.206004728741926]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-357.46364178007343]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 2.3399671478555746e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0008847975558072403
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
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
        answer = 'Forsterite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg2.000SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1984446.1221989095]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 33.450533005483145
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 357.46364178007343
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
        answer = 'Forsterite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 4.206004728741926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestMgPerovskite:
    phase = db.MgPerovskite()

    def test_name(self):
        test = self.phase.name()
        answer = 'MgPerovskite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 2455273.6800218695
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.302105774495801
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 100.3887
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'MgPv'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 3.207220234066898e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 4.0728657181349043e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 100.3887
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 133.48390635877956
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 122.9723459619135
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[2.4482783601974316]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-245.96599745551111]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -9.971509001699658e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.07851994491692915
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 7.852167895453328e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-9.971509001699658e-07]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.07851994491692915]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[7.852167895453328e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 2.1533239202025615e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 1.544309042514691e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 3.990400107339729e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -1.7460043586313647e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.010683143092153763
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1363449.0604146551]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 2.4482783601974316
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -245.96599745551111
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[2.4482783601974316]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-245.96599745551111]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 1.670029443405111e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.00131508195777913
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
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
        answer = 'MgPerovskite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg1.000SiO3'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1363449.0604146551]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 41.00379337254141
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 245.96599745551111
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
        answer = 'MgPerovskite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 2.4482783601974316
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestMgRingwoodite:
    phase = db.MgRingwoodite()

    def test_name(self):
        test = self.phase.name()
        answer = 'MgRingwoodite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1916858.2072982239
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.5245310300166635
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'MgRi'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 2.4539534710325686e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 5.216869960399844e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 179.7983142829723
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 172.1994273285761
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[3.8723840789560926]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-345.72095455231477]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -2.0201724176636655e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10576371428410136
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 9.50265035172556e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-2.0201724176636655e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10576371428410136]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[9.50265035172556e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 5.822290435919738e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.0705894836176042e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.666994947427094e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -3.4060101483156335e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.009424800177840767
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1966692.8925920525]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 3.8723840789560926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -345.72095455231477
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[3.8723840789560926]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-345.72095455231477]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 1.8954153102070262e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0008915807782899598
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
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
        answer = 'MgRingwoodite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg2.000SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1966692.8925920525]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 36.33242393609047
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 345.72095455231477
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
        answer = 'MgRingwoodite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 3.8723840789560926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestMgWadsleyite:
    phase = db.MgWadsleyite()

    def test_name(self):
        test = self.phase.name()
        answer = 'MgWadsleyite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1757071.3405069287
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.6238090125529325
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'MgWa'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 2.8670619149850206e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 5.6912885489982e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 140.6931
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 182.13818429275685
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 172.37672081299118
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[3.975598095355361]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-352.3823880256589]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -2.262627591551502e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10714010840750403
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.00011398285888480343
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-2.262627591551502e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10714010840750403]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.00011398285888480343]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 7.2419287413493e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.64020860816539e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.645682346622832e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -4.2741988691355865e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.011163508514915876
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1975655.3236001502]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 3.975598095355361
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -352.3823880256589
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[3.975598095355361]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-352.3823880256589]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 2.014099538604103e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0010146292935959437
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
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
        answer = 'MgWadsleyite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg2.000SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-1975655.3236001502]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 35.38916576209499
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 352.3823880256589
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
        answer = 'MgWadsleyite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 3.975598095355361
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
        test = self.phase.K(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 1429347.697586451
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 4.405206893378452
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([0.9, 0.1]),)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Ol'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.6335852373441714e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 6.996198347599865e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([0.86137754, 0.13862246]),)
        answer =[0.9, 0.09999999999999999]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([0.9, 0.1]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_mole_frac(self):
        test = self.phase.conv_moles_to_mole_frac(array([0.9, 0.1]),)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_tot_moles(self):
        test = self.phase.conv_moles_to_tot_moles(array([0.9, 0.1]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 179.72393146424304
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 172.59103253316195
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[2989.0192113025114, -26901.172901721206, 242110.55611548852]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[4.206004728741926, 4.469742912034868]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-359.2156739178061, -442.1105522473071]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -2.9610559797436797e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -0.10571995968484885
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.00011146329660418941
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[1.8476583595896083, -16.62892523630648, 149.66032712675832]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-6175.040557191563, 25685.172901723534, 37845.172901724, -2761712.1172703765]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-2.942229018403993e-06, -3.1304986318008575e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-0.10567913485286703, -0.10608738317268504]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[0.00011125271765010549, 0.00011335850719094473]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 1.1197499545013253e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.4357969646912492e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.7176569415573865e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -5.145848806447156e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.008519791678373273
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1987348.5768330551, -1394215.519369004]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 4.232378547071221
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -367.5051617507562
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[4.206004728741926, 4.469742912034868]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-359.2156739178061, -442.1105522473071]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[[3123.0409138064733, -19406.070497187313], [-28107.368224257923, 174654.63447468376]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 2.42996140339664e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = -0.0009147125569269987
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-1.3312022689960887, 8.271875326368122]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-8.661172052328183, 53.819120554167895]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-0.02755636811708456, 0.17123080903695057]
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
        answer = 'Forsterite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 2
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 'Mg1.800Fe0.200SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -1928035.2710866502
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
        test = self.phase.mu(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1987348.5768330551, -1394215.519369004]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 34.73259737168928
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 367.5051617507562
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
        answer = 'Forsterite_stixrude'
        assert(test == answer)

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 2
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([0.9, 0.1]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 4.232378547071221
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([0.9, 0.1]),)
        answer =[0.8613775369639085, 0.1386224630360915]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestPericlase:
    phase = db.Periclase()

    def test_name(self):
        test = self.phase.name()
        answer = 'Periclase'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1629319.3004791553
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.1904091850190825
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 40.3044
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Pe'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 3.636368711886252e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 6.137532402064573e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 40.3044
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 53.43230243192754
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 49.365418571070265
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[1.110377932705695]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-105.63270649432897]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -6.814980540518679e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.031430766136427965
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 4.0377435728599276e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-6.814980540518679e-07]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.031430766136427965]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[4.0377435728599276e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 2.1710009562172374e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 1.0148203879193656e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 1.6002932955665014e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -1.4575289337765742e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.004225780111797441
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-563891.0645437283]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 1.110377932705695
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -105.63270649432897
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[1.110377932705695]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-105.63270649432897]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 2.2277960832938898e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0013199258990514748
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'MgO'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 40.3044
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'Periclase_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Mg1.000O'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-563891.0645437283]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 36.29791156042602
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 105.63270649432897
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'MgO'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 40.3044
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'Periclase_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 1.110377932705695
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestRingwoodite:
    phase = db.Ringwoodite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Ringwoodite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 1940137.2859308969
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 4.536620086155425
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([0.9, 0.1]),)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Ri'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.447090647664383e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.154274428163419e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([0.86137754, 0.13862246]),)
        answer =[0.9, 0.09999999999999999]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([0.9, 0.1]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_mole_frac(self):
        test = self.phase.conv_moles_to_mole_frac(array([0.9, 0.1]),)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_tot_moles(self):
        test = self.phase.conv_moles_to_tot_moles(array([0.9, 0.1]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 180.00218729690516
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 172.3033428528836
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[2959.0192113020457, -26631.17290172167, 239680.55611548852]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[3.8723840789560926, 4.1287418291048485]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-347.47298669004743, -429.37622967158194]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -2.0091464053795864e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -0.1058836395864148
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 9.53880792906244e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[1.8476583595896083, -16.62892523630648, 149.66032712675832]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-6085.040557195288, 25175.17290171981, 39735.172901724465, -2754422.1172703775]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-2.0201724176636655e-06, -1.909912294822873e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-0.10576371428410136, -0.10696296730723558]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[9.50265035172556e-05, 9.864226125094362e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.733553199929479e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.0677865812993057e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.67653674074258e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -3.3813794311073833e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.009382514993790933
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1969580.3472261978, -1405190.4933938738]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 3.8980198539709683
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -355.66331098820086
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[3.8723840789560926, 4.1287418291048485]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-347.47298669004743, -429.37622967158194]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[[3091.6957732132246, -19211.296869532394], [-27825.26195891873, 172901.67182578973]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 1.9437717116289162e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = -0.0009228429031118908
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-1.4317874155442143, 8.896895138389628]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-8.55756222099734, 53.175305840720284]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-0.026785232401939026, 0.16643909657975844]
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
        answer = 'MgRingwoodite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 2
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 'Mg1.800Fe0.200SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -1913141.3618429657
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
        test = self.phase.mu(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1969580.3472261978, -1405190.4933938738]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 37.711839730689796
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 355.66331098820086
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
        answer = 'MgRingwoodite_stixrude'
        assert(test == answer)

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 2
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([0.9, 0.1]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 3.8980198539709683
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([0.9, 0.1]),)
        answer =[0.8613775369639085, 0.1386224630360915]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestWadsleyite:
    phase = db.Wadsleyite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Wadsleyite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 1758489.5350845163
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 4.6204077424780685
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([0.9, 0.1]),)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Wa'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.849871956496243e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.68669861292031e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([0.86137754, 0.13862246]),)
        answer =[0.9, 0.09999999999999999]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([0.9, 0.1]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_mole_frac(self):
        test = self.phase.conv_moles_to_mole_frac(array([0.9, 0.1]),)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_tot_moles(self):
        test = self.phase.conv_moles_to_tot_moles(array([0.9, 0.1]),)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 182.16806499035764
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 172.46224757114186
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[2811.0192113025114, -25299.172901721206, 227692.55611548852]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[3.975598095355361, 4.195007402310355]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-354.13442016339155, -431.41792482994674]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -2.273279963453824e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -0.10715768528844567
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.00011392474365341653
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[1.8476583595896083, -16.62892523630648, 149.66032712675832]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-5641.040557193425, 22659.172901716083, 49059.17290172214, -2718458.1172703765]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-2.262627591551502e-06, -2.369151310574721e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-0.10714010840750403, -0.1073158772169204]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[0.00011398285888480343, 0.0001134017065709345]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 7.26575851177986e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.6164319774048037e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.656943163009776e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -4.26560116378943e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.010989651517279486
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1978468.7782342958, -1396605.6743665363]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 3.9975390260508608
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -361.86277063004707
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[3.975598095355361, 4.195007402310355]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-354.13442016339155, -431.41792482994674]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[[2937.059746286531, -18250.413639768107], [-26433.53771657844, 164253.72275791084]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 2.0911696438722132e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = -0.0010479833959913729
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-1.4379460625871288, 8.935163973790942]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-8.07487439770721, 50.17596187263706]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-0.022924718579904937, 0.14245048885625583]
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
        answer = 'MgWadsleyite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 2
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 'Mg1.800Fe0.200SiO4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -1920282.4678475198
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
        test = self.phase.mu(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1978468.7782342958, -1396605.6743665363]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 36.772999348357004
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 361.86277063004707
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
        answer = 'MgWadsleyite_stixrude'
        assert(test == answer)

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 2
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([0.9, 0.1]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 3.9975390260508608
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([0.9, 0.1]),)
        answer =[0.8613775369639085, 0.1386224630360915]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestWuestite:
    phase = db.Wuestite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Wuestite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1805805.4394082627
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 5.32196804073994
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([1.]),)
        answer = 71.8464
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Wu'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, array([1.]))
        answer = 3.415191676976622e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 5.537695136900718e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 71.8464
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],)
        answer = 1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_moles_to_elm(self):
        test = self.phase.conv_moles_to_elm(array([1.]),)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 54.05994713620336
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 49.704773175038255
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[1.216341324578678]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-131.82994540737673]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -6.735727437930723e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.031799968903649035
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 4.1540387680638215e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([1.]))
        answer =[-6.735727437930723e-07]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.031799968903649035]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[4.1540387680638215e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 2.358119688003987e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 1.110575181748639e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 1.6305900476593658e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -1.613064198088169e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.004079938093439817
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-261701.67654927325]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 1.216341324578678
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -131.82994540737673
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[1.216341324578678]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-131.82994540737673]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 3.270985305227852e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0020172728028107186
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_elements(self):
        test = self.phase.endmember_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_formula(self):
        test = self.phase.endmember_formula(0,)
        answer = 'FeO'
        assert(test == answer)

    def test_endmember_mw(self):
        test = self.phase.endmember_mw(0,)
        answer = 71.8464
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_name(self):
        test = self.phase.endmember_name(0,)
        answer = 'Wuestite_stixrude'
        assert(test == answer)

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 1
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 100000.0, array([1.]))
        answer = 'Fe1.000O'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 100000.0, array([1.]))
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

    def test_mu(self):
        test = self.phase.mu(1700.0, 100000.0, array([1.]))
        answer =[-261701.67654927325]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 59.06763056404952
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 131.82994540737673
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_elements(self):
        test = self.phase.species_elements(0,)
        answer =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_formula(self):
        test = self.phase.species_formula(0,)
        answer = 'FeO'
        assert(test == answer)

    def test_species_mw(self):
        test = self.phase.species_mw(0,)
        answer = 71.8464
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_name(self):
        test = self.phase.species_name(0,)
        answer = 'Wuestite_stixrude'
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
        test = self.phase.v(1700.0, 100000.0, array([1.]))
        answer = 1.216341324578678
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([1.]),)
        answer =[1.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

