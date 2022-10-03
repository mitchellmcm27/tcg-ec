import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_fa_binary as db
import pytest
class TestLiquid:
    phase = db.Liquid()

    def test_name(self):
        test = self.phase.name()
        answer = 'Liquid'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 292095.07071152935
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 11.404159120571604
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([0.68479816, 0.31520184]),)
        answer = 160.57729297347703
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Liq'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 0.00010455569797812429
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 3.423542881309324e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([0.6, 0.4]),)
        answer =[0.6847981584319788, 0.3152018415680212]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 261.291783279705
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 233.52355780418046
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[6505.928242148118, -14134.586450860506, 30708.382678208596]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[4.980615135732, 5.408228374505001]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[-414.3607676410934, -527.3917927989655]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = -1.7512789940656142e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = -0.15370104898806175
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 0.0005348441772954491
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[3.827016613028305, -8.31446261815324, 18.063754516593296]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[-16006.43311652616, 14134.586450860506, 14134.586450860506, -128132.88542046706]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[-1.3260939e-05, -2.6750235e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[-0.15941176470588236, -0.14129411764705882]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[0.0005175065000000001, 0.0005725115]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 7.437011262801502e-10
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 9.041238175768337e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = -8.076009207840106e-09
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 2.7755575615628914e-17
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[-2549061.2183886603, -2017782.6286599443]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 5.115399616072116
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = -449.988354925176
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[4.980615135732, 5.408228374505001]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[-414.3607676410934, -527.3917927989655]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 1000.0, array([0.6, 0.4]))
        answer = 0.000107468289775842
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 1000.0, array([0.6, 0.4]))
        answer = -0.0032821034926634658
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 2
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 'Mg1.370Fe0.630Si2O4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = -2381601.228520508
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
        test = self.phase.mu(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer =[-2549061.2183886603, -2017782.6286599443]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, array([0.6, 0.4]))
        answer = 31.390957701321696
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 449.988354925176
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 2
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([0.68479816, 0.31520184]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0, array([0.68479816, 0.31520184]))
        answer = 5.115399616072116
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([0.68479816, 0.31520184]),)
        answer =[0.6, 0.4000000000000001]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))


class TestOlivine:
    phase = db.Olivine()

    def test_name(self):
        test = self.phase.name()
        answer = 'Olivine'
        assert(test == answer)

    def test_K(self):
        test = self.phase.K(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 1344309.2874781073
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Kp(self):
        test = self.phase.Kp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = -1.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mass(self):
        test = self.phase.Mass(array([0.9287517, 0.0712483]),)
        answer = 145.1877280197317
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_abbrev(self):
        test = self.phase.abbrev()
        answer = 'Ol'
        assert(test == answer)

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 5.104822851920045e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 7.438764347719241e-07
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_c_to_x(self):
        test = self.phase.c_to_x(array([0.9, 0.1]),)
        answer =[0.9287516958383784, 0.07124830416162155]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_cp(self):
        test = self.phase.cp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 190.8149112047725
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_cv(self):
        test = self.phase.cv(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 163.19882862636067
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdn2(self):
        test = self.phase.d2gdn2(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[1084.3213736913574, -14134.58645086051, 184250.29606925166]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndp(self):
        test = self.phase.d2gdndp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[4.619242458468521, 4.870800598538133]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdndt(self):
        test = self.phase.d2gdndt(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[-373.6307385824196, -474.7437423316388]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdp2(self):
        test = self.phase.d2gdp2(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = -3.4494781762950887e-06
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdt2(self):
        test = self.phase.d2gdt2(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = -0.11224406541457206
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d2gdtdp(self):
        test = self.phase.d2gdtdp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 0.00023671908664440544
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dp(self):
        test = self.phase.d3gdn2dp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[0.0, 0.0, 0.0]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn2dt(self):
        test = self.phase.d3gdn2dt(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[0.6378361021713868, -8.31446261815324, 108.38252709955978]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdn3(self):
        test = self.phase.d3gdn3(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[-2251.8254316113266, 14134.586450860505, 14134.586450860503, -2770280.9144442994]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndp2(self):
        test = self.phase.d3gdndp2(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[-3.4548158e-06, -3.3799e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndt2(self):
        test = self.phase.d3gdndt2(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[-0.11173412815267972, -0.11889131270454778]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdndtdp(self):
        test = self.phase.d3gdndtdp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[0.00023758300657999998, 0.0002254575259]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdp3(self):
        test = self.phase.d3gdp3(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt2dp(self):
        test = self.phase.d3gdt2dp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 7.738985073437966e-08
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdt3(self):
        test = self.phase.d3gdt3(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 5.767210801118094e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_d3gdtdp2(self):
        test = self.phase.d3gdtdp2(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dcpdt(self):
        test = self.phase.dcpdt(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 0.014201481795564463
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdn(self):
        test = self.phase.dgdn(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[-2567390.214882203, -2025944.878949086]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdp(self):
        test = self.phase.dgdp(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 4.637165549346531
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dgdt(self):
        test = self.phase.dgdt(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = -380.83486862823906
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dP(self):
        test = self.phase.dmu_dP(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[4.619242458468521, 4.870800598538133]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dmu_dT(self):
        test = self.phase.dmu_dT(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[-373.6307385824196, -474.7437423316388]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dP(self):
        test = self.phase.drho_dP(1700.0, 1000.0, array([0.9, 0.1]))
        answer = 2.329046231855436e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_drho_dT(self):
        test = self.phase.drho_dT(1700.0, 1000.0, array([0.9, 0.1]))
        answer = -0.001598298839940969
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_endmember_number(self):
        test = self.phase.endmember_number()
        answer = 2
        assert(test == answer)

    def test_formula(self):
        test = self.phase.formula(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 'Mg1.858Fe0.142Si2O4'
        assert(test == answer)

    def test_g(self):
        test = self.phase.g(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = -2528813.1529007484
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
        test = self.phase.mu(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer =[-2567390.214882203, -2025944.878949086]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, array([0.9, 0.1]))
        answer = 31.309584804491514
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 380.83486862823906
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_species_number(self):
        test = self.phase.species_number()
        answer = 2
        assert(test == answer)

    def test_test_moles(self):
        test = self.phase.test_moles(array([0.9287517, 0.0712483]),)
        answer = 1
        assert(test == answer)

    def test_v(self):
        test = self.phase.v(1700.0, 1000.0, array([0.9287517, 0.0712483]))
        answer = 4.637165549346531
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_x_to_c(self):
        test = self.phase.x_to_c(array([0.9287517, 0.0712483]),)
        answer =[0.9, 0.1]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))
