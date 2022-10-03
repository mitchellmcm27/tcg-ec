import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_fa_binary as db
import pytest
class Testfo_fa_binary:
    phase = db.fo_fa_binary()

    def test_name(self):
        test = self.phase.name()
        answer = 'fo_fa_binary'
        assert(test == answer)

    def test_A(self):
        test = self.phase.A(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[-18328.9964935421, -8162.2502891412005]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_C_to_X(self):
        test = self.phase.C_to_X([[0.9, 0.1], [0.6, 0.4]],)
        answer =[[0.9287516958383784, 0.07124830416162155], [0.6847981584319788, 0.3152018415680212]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Cp(self):
        test = self.phase.Cp(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[1.314263359633535, 1.6272025666970438]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_i(self):
        test = self.phase.Gamma_i(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[0.8062166068802055, -0.8062166068802057]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_ik(self):
        test = self.phase.Gamma_ik(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[[0.5578122268749602, 0.24840438000524526], [-0.5578122268749603, -0.24840438000524534]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_M(self):
        test = self.phase.M()
        answer =[[140.69310000000002, 203.7731], [140.69310000000002, 203.7731]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mu(self):
        test = self.phase.Mu(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[[-2567390.214882203, -2025944.878949086], [-2549061.2183886603, -2017782.6286599443]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_X_to_C(self):
        test = self.phase.X_to_C([[0.9287516958383784, 0.07124830416162155], [0.6847981584319788, 0.3152018415680212]],)
        answer =[[0.9, 0.1], [0.6, 0.4000000000000001]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[5.104822851920045e-05, 0.00010455569797812429]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dC(self):
        test = self.phase.dA_dC(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[[[11189.61460721827, -11189.614607218267], [-18563.531996306803, 18563.531996306803]], [[-145861.3459578985, 145861.3459578985], [40330.57821561199, -40330.57821561199]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dP(self):
        test = self.phase.dA_dP(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[-0.3613726772634802, -0.5374277759668686]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dT(self):
        test = self.phase.dA_dT(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[40.73002905867389, 52.64805046732681]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dAj_dCik(self):
        test = self.phase.dAj_dCik(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], 0, 0, 1)
        answer = -11189.614607218267
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dC(self):
        test = self.phase.dGamma_i_dC(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0, 0, 1)
        answer = -4.098507978235309
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dP(self):
        test = self.phase.dGamma_i_dP(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0)
        answer = 2.735348236382119e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dPhi(self):
        test = self.phase.dGamma_i_dPhi(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0, 0)
        answer = -0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dT(self):
        test = self.phase.dGamma_i_dT(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0)
        answer = -0.0030370824116391585
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dC(self):
        test = self.phase.dGamma_ik_dC(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0, 1, 0, 1)
        answer = -4.4390451071569155
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dP(self):
        test = self.phase.dGamma_ik_dP(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0, 1)
        answer = 1.6355711814455288e-05
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dPhi(self):
        test = self.phase.dGamma_ik_dPhi(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0, 1, 0)
        answer = -0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dT(self):
        test = self.phase.dGamma_ik_dT(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 0, 1)
        answer = -0.001662422150595961
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_get_parameter(self):
        test = self.phase.get_parameter('R',)
        answer = None
        assert(test == answer)

    def test_list_parameters(self):
        test = self.phase.list_parameters()
        answer = None
        assert(test == answer)

    def test_nu(self):
        test = self.phase.nu()
        answer =[[[-1.0, 0.0], [1.0000000000000002, 0.0]], [[0.0, -1.0], [0.0, 1.0000000000000002]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_nu_m(self):
        test = self.phase.nu_m()
        answer =[[[-0.9999999999999998, 0.0], [1.0, 0.0]], [[0.0, -0.9999999999999997], [0.0, 1.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_report(self):
        test = self.phase.report()
        answer = None
        assert(test == answer)

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, 0)
        answer = 28.248137261326928
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_rho_phases(self):
        test = self.phase.rho_phases(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[31.309584804491514, 31.390957701321696]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[2.6230513682016, 2.8023162340861094]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_set_parameter(self):
        test = self.phase.set_parameter('R', 8.314)
        answer = None
        assert(test == answer)

    def test_zero_C(self):
        test = self.phase.zero_C()
        answer =[[0.0, 0.0], [0.0, 0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

