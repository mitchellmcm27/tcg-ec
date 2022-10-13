import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_h2o_hydration as db
import pytest
class Testfo_h2o_hydration:
    phase = db.fo_h2o_hydration()

    def test_name(self):
        test = self.phase.name()
        answer = 'fo_h2o_hydration'
        assert(test == answer)

    def test_A(self):
        test = self.phase.A(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[48123.598109198036]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_C_to_X(self):
        test = self.phase.C_to_X([[1.0], [1.0], [1.0], [1.0]],)
        answer =[[1.0], [1.0], [1.0], [1.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Cp(self):
        test = self.phase.Cp(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[0.993790544249227, 1.1941436828632443, 1.6059992821104447, 4.308570879729805]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_i(self):
        test = self.phase.Gamma_i(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[-1.6941351875926878, 1.668403273325389, 0.3511240326132275, -0.32539211834592907]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_ik(self):
        test = self.phase.Gamma_ik(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[[-1.6941351875926878], [1.668403273325389], [0.3511240326132275], [-0.32539211834592907]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_M(self):
        test = self.phase.M()
        answer =[[140.69332], [277.112709], [58.31979199999999], [18.015287]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mu(self):
        test = self.phase.Mu(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[[-2217084.332036725], [-4463217.900778551], [-954530.8351882802], [-247654.0271157971]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_X_to_C(self):
        test = self.phase.X_to_C([[1.0], [1.0], [1.0], [1.0]],)
        answer =[[1.0], [1.0], [1.0], [1.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[3.156827010221178e-05, 2.8726266313116093e-05, 3.299333900718681e-05, 0.0010244809642597453]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dC(self):
        test = self.phase.dA_dC(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[[[0.0], [0.0], [0.0], [0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dP(self):
        test = self.phase.dA_dP(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[0.2843321061990641]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dT(self):
        test = self.phase.dA_dT(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[-29.781503411911192]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dAj_dCik(self):
        test = self.phase.dAj_dCik(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], 0, 1, 0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dC(self):
        test = self.phase.dGamma_i_dC(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[[[-0.0], [-0.0], [-0.0], [-0.0]], [[0.0], [0.0], [0.0], [0.0]], [[0.0], [0.0], [0.0], [0.0]], [[-0.0], [-0.0], [-0.0], [-0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dP(self):
        test = self.phase.dGamma_i_dP(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[-1.0009580434554145e-05, 9.857546719960464e-06, 2.074571304986147e-06, -1.9225375903924685e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dPhi(self):
        test = self.phase.dGamma_i_dPhi(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[[-0.941186215329271, -0.0, -0.0, -8.470675937963438], [0.9268907074029941, 0.0, 0.0, 8.342016366626945], [0.19506890700734864, 0.0, 0.0, 1.7556201630661377], [-0.1807733990810717, -0.0, -0.0, -1.6269605917296452]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dT(self):
        test = self.phase.dGamma_i_dT(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[0.0022136051783073893, -0.0021799831279025496, -0.0004587886388956129, 0.00042516658849077335]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dC(self):
        test = self.phase.dGamma_ik_dC(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1], 1)
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dP(self):
        test = self.phase.dGamma_ik_dP(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[[-1.0009580434554145e-05], [9.857546719960464e-06], [2.074571304986147e-06], [-1.9225375903924685e-06]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dPhi(self):
        test = self.phase.dGamma_ik_dPhi(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1], 1)
        answer =[[0.9268907074029941, 0.0, 0.0, 8.342016366626945]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dT(self):
        test = self.phase.dGamma_ik_dT(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]], [0.9, 0.0, 0.0, 0.1])
        answer =[[0.0022136051783073893], [-0.0021799831279025496], [-0.0004587886388956129], [0.00042516658849077335]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dMu_dC(self):
        test = self.phase.dMu_dC(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[[[0.0]], [[0.0]], [[0.0]], [[0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dC(self):
        test = self.phase.ds_dC(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[[0.0], [0.0], [0.0], [0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

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
        answer =[[[-0.4], [0.2], [0.2], [-0.6]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_nu_m(self):
        test = self.phase.nu_m()
        answer =[[[-0.8388770890152949], [0.8261355359837358], [0.17386446401626415], [-0.1611229109847051]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_report(self):
        test = self.phase.report()
        answer = None
        assert(test == answer)

    def test_rho(self):
        test = self.phase.rho(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[32.102246729760815, 25.76031332955861, 23.533685517229163, 9.173042157504625]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(423.15, 10.0, [[1.0], [1.0], [1.0], [1.0]])
        answer =[0.9907161200838106, 1.1772799456824503, 1.5972083523946317, 5.356949965856523]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_set_parameter(self):
        test = self.phase.set_parameter('R', 8.314)
        answer = None
        assert(test == answer)

    def test_zero_C(self):
        test = self.phase.zero_C()
        answer =[[0.0], [0.0], [0.0], [0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

