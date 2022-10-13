import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_MgFeSiO4_stixrude as db
import pytest
class TestFayalite:
    phase = db.Fayalite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Fayalite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1427805.4194400304
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.692110212623922
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
        answer = 2.5361303641362687e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 7.003755458444674e-07
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
        answer = 203.77710000000002
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
        answer = 180.3485513935645
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 173.3703477497621
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[4.469742912034869]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-403.8210368856753]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -3.130498631800861e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.106087383172685
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.00011335850719094499
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
        answer =[-3.130498631800861e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.106087383172685]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.00011335850719094499]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 1.248009216806823e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 3.253008732157926e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.76184593765122e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -6.620961655580267e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.00813600223261425
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1335279.3432542302]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.469742912034869
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -403.8210368856753
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.469742912034869]
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
        answer = 3.193035940810487e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0011562304611169574
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
        answer = 45.59034020756009
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
        answer = 4.469742912034869
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
        answer = 2161744.20170836
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.585612787379428
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
        answer = 2.3891603140594054e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 4.625894216391239e-07
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
        answer = 203.77710000000002
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
        answer = 181.8370444223006
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 173.1761741075603
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[4.128741829104848]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-391.08671430994997]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -1.90991229482287e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10696296730723565
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 9.864226125094341e-05
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
        answer =[-1.90991229482287e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10696296730723565]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[9.864226125094341e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 4.9349180760171315e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.0425604604346447e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.762412880581907e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -3.1597029762331156e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.009001948337343232
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1347469.3172790997]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.128741829104848
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -391.08671430994997
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.128741829104848]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-391.08671430994997]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 2.283144229745543e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0011791877050827133
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
        answer = -1347469.3172790997
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
        answer =[-1347469.3172790997]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 49.35573800316328
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 391.08671430994997
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
        answer = 4.128741829104848
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
        answer = 1770679.3920615825
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.590644530028838
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
        answer = 2.7032540278350895e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 5.647549773738039e-07
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
        answer = 203.77710000000002
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
        answer = 182.43699126876507
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 173.20925224979237
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[4.195007402310356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-393.128409468315]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -2.3691513105747253e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10731587721692062
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.00011340170657093485
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
        answer =[-2.3691513105747253e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10731587721692062]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.00011340170657093485]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 7.480226445654928e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.4024423005595503e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.7582905104922704e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -4.188221815674043e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.00942493853855203
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1344878.4982517632]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.195007402310356
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -393.128409468315
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.195007402310356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-393.128409468315]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 2.7433594380886672e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0013131353857735095
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
        answer = -1344878.4982517632
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
        answer =[-1344878.4982517632]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 48.576100220412464
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 393.128409468315
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
        answer = 4.195007402310356
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
        answer = 4.371247101982013
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
        answer = 2.645092548038645e-05
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 179.6545292498743
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 172.50308587748793
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
        answer =[-357.4636417800734]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -2.942229018403993e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10567913485286723
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.00011125271765010572
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
        answer =[-0.10567913485286723]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.00011125271765010572]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 1.10549892535627e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.3449956571950204e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.712747053102467e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -4.981947378765714e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.008562434950125292
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1984446.1221989086]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 4.206004728741926
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -357.4636417800734
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[4.206004728741926]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-357.4636417800734]
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
        answer = -0.000884797555807242
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
        answer = -1984446.1221989086
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
        answer =[-1984446.1221989086]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 33.450533005483145
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 357.4636417800734
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


class TestMgRingwoodite:
    phase = db.MgRingwoodite()

    def test_name(self):
        test = self.phase.name()
        answer = 'MgRingwoodite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([1.]))
        answer = 1916858.2072982246
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.524531030016663
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
        answer = 2.4539534710325676e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 5.216869960399842e-07
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 179.7983142829724
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 172.1994273285762
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([1.]))
        answer =[0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([1.]))
        answer =[3.872384078956092]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([1.]))
        answer =[-345.7209545523147]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([1.]))
        answer = -2.0201724176636647e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.1057637142841014
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 9.502650351725554e-05
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
        answer =[-2.0201724176636647e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.1057637142841014]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[9.502650351725554e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 5.822290435919733e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.070589483617603e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.666994947427094e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -3.406010148315631e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.009424800177840809
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([1.]))
        answer =[-1966692.892592052]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([1.]))
        answer = 3.872384078956092
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([1.]))
        answer = -345.7209545523147
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([1.]))
        answer =[3.872384078956092]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([1.]))
        answer =[-345.7209545523147]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([1.]))
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([1.]))
        answer = 1.895415310207026e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0008915807782899595
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
        answer = -1966692.892592052
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
        answer =[-1966692.892592052]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([1.]))
        answer = 36.33242393609047
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([1.]))
        answer = 345.7209545523147
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
        answer = 3.872384078956092
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
        answer = 1757071.3405069292
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([1.]))
        answer = 4.62380901255293
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
        answer = 2.8670619149850216e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([1.]))
        answer = 5.691288548998199e-07
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
        test = self.phase.cp(1700.0, 100000.0, array([1.]))
        answer = 182.1381842927569
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([1.]))
        answer = 172.37672081299124
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
        answer = -2.2626275915515014e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([1.]))
        answer = -0.10714010840750406
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([1.]))
        answer = 0.00011398285888480347
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
        answer =[-2.2626275915515014e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([1.]))
        answer =[-0.10714010840750406]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([1.]))
        answer =[0.00011398285888480347]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([1.]))
        answer = 7.241928741349294e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([1.]))
        answer = 2.640208608165391e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([1.]))
        answer = 5.645682346622833e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([1.]))
        answer = -4.2741988691355834e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([1.]))
        answer = 0.011163508514915904
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
        answer = 2.0140995386041028e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([1.]))
        answer = -0.0010146292935959442
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
        answer = 4.405206893378453
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
        answer = 2.633585237344177e-05
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
        test = self.phase.conv_elm_to_moles(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
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
        answer = 179.7239314642433
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 172.59103253316218
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[2989.0192113025114, -26901.172901721206, 242110.55611548852]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[4.206004728741926, 4.469742912034869]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-359.21567391780604, -442.1105522473071]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -2.9610559797436797e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -0.105719959684849
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.00011146329660418965
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
        answer =[-6175.040557194501, 25685.17290172167, 37845.17290172167, -2761712.1172703747]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-2.942229018403993e-06, -3.130498631800861e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-0.10567913485286723, -0.106087383172685]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[0.00011125271765010572, 0.00011335850719094499]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 1.1197499545013254e-11
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.435796964691311e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.7176569415573425e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -5.145848806447169e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.008519791678374175
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1987348.5768330542, -1394215.519369004]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 4.232378547071221
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -367.50516175075614
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[4.206004728741926, 4.469742912034869]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-359.21567391780604, -442.1105522473071]
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
        answer = -0.0009147125569270006
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-1.3312022689960887, 8.271875326368114]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-8.661172052328242, 53.819120554167895]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-0.02755636811708456, 0.1712308090369512]
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
        answer = -1928035.2710866493
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
        answer =[-1987348.5768330542, -1394215.519369004]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 34.73259737168928
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 367.50516175075614
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


class TestRingwoodite:
    phase = db.Ringwoodite()

    def test_name(self):
        test = self.phase.name()
        answer = 'Ringwoodite'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 1940137.2859308978
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 4.536620086155424
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
        answer = 2.4470906476643816e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.154274428163415e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([0.86137754, 0.13862246]),)
        answer =[0.9, 0.09999999999999999]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_moles(self):
        test = self.phase.conv_elm_to_moles(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
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
        answer = 180.00218729690522
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 172.30334285288365
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[2959.0192113020457, -26631.172901721206, 239680.556115489]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[3.872384078956092, 4.128741829104848]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-347.4729866900474, -429.37622967158177]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -2.009146405379585e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -0.10588363958641483
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 9.538807929062433e-05
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
        answer =[-6085.0405571926385, 25175.17290171981, 39735.17290172167, -2754422.1172703784]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-2.0201724176636647e-06, -1.90991229482287e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-0.1057637142841014, -0.10696296730723565]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[9.502650351725554e-05, 9.864226125094341e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.7335531999294725e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.0677865812993074e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.676536740742575e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -3.381379431107379e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.009382514993791058
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1969580.3472261974, -1405190.4933938738]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 3.898019853970968
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -355.6633109882008
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[3.872384078956092, 4.128741829104848]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-347.4729866900474, -429.37622967158177]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[[3091.6957732131755, -19211.29686953209], [-27825.261958918727, 172901.67182578973]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 1.943771711628916e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = -0.0009228429031118907
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-1.431787415544219, 8.89689513838963]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-8.55756222099734, 53.175305840720206]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-0.026785232401938565, 0.16643909657975844]
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
        answer = -1913141.3618429652
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
        answer =[-1969580.3472261974, -1405190.4933938738]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 37.7118397306898
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 355.6633109882008
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
        answer = 3.898019853970968
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
        answer = 4.620407742478067
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
        answer = 2.8498719564962452e-05
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
        test = self.phase.conv_elm_to_moles(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_grams(self):
        test = self.phase.conv_elm_to_tot_grams(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
        answer = 147.0015
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_conv_elm_to_tot_moles(self):
        test = self.phase.conv_elm_to_tot_moles(array([0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 4. , 0. , 0. , 0. , 1.8,
       0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. ]),)
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
        answer = 182.16806499035772
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 172.46224757114192
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[2811.019211302977, -25299.17290172074, 227692.55611548852]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[3.975598095355361, 4.195007402310356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-354.13442016339155, -431.4179248299468]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -2.273279963453824e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -0.10715768528844571
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.00011392474365341661
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
        answer =[-5641.040557196364, 22659.172901717946, 49059.17290171981, -2718458.1172703803]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-2.2626275915515014e-06, -2.3691513105747253e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-0.10714010840750406, -0.10731587721692062]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[0.00011398285888480347, 0.00011340170657093485]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 7.2657585117798576e-12
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 2.616431977404807e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 5.6569431630097766e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 100000.0, array([0.9, 0.1]))
        answer = -4.2656011637894293e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 100000.0, array([0.9, 0.1]))
        answer = 0.010989651517279514
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-1978468.778234296, -1396605.674366537]
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
        answer =[3.975598095355361, 4.195007402310356]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 100000.0, array([0.9, 0.1]))
        answer =[-354.13442016339155, -431.4179248299468]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dc(self):
        test = self.phase.dmu_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[[2937.0597462865317, -18250.413639768107], [-26433.53771657839, 164253.72275791055]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = 2.091169643872213e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer = -0.0010479833959913733
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dc(self):
        test = self.phase.drho_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-1.4379460625871243, 8.935163973790937]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dc(self):
        test = self.phase.ds_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-8.07487439770721, 50.1759618726371]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dv_dc(self):
        test = self.phase.dv_dc(1700.0, 100000.0, array([0.86137754, 0.13862246]))
        answer =[-0.022924718579905398, 0.14245048885625614]
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
        answer = -1920282.46784752
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
        answer =[-1978468.778234296, -1396605.674366537]
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

