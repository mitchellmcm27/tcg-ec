import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_Mg2SiO4_stixrude as db
import pytest
class TestMg2SiO4_stixrude:
    phase = db.Mg2SiO4_stixrude()

    def test_name(self):
        test = self.phase.name()
        answer = 'Mg2SiO4_stixrude'
        assert(test == answer)

    def test_A(self):
        test = self.phase.A(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[-8790.79859875841, -8962.431008098181]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_C_to_X(self):
        test = self.phase.C_to_X([[1.0], [1.0], [1.0]],)
        answer =[[1.0], [1.0], [1.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Cp(self):
        test = self.phase.Cp(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[1.2769249469225876, 1.294577945135596, 1.2779469233599403]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_i(self):
        test = self.phase.Gamma_i(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[0.06329197443935008, 0.0012357186827634814, -0.06452769312211357]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_ik(self):
        test = self.phase.Gamma_ik(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[[0.06329197443935008], [0.0012357186827634814], [-0.06452769312211357]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_M(self):
        test = self.phase.M()
        answer =[[140.69332], [140.69332], [140.69332]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mu(self):
        test = self.phase.Mu(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[[-1984446.1221989086], [-1975655.3236001502], [-1966692.892592052]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_X_to_C(self):
        test = self.phase.X_to_C([[1.0], [1.0], [1.0]],)
        answer =[[1.0], [1.0], [1.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[2.645092548038645e-05, 2.8670619149850216e-05, 2.4539534710325676e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dC(self):
        test = self.phase.dA_dC(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[[[0.0], [0.0], [0.0]], [[0.0], [0.0], [0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dP(self):
        test = self.phase.dA_dP(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[0.2304066333865653, 0.10321401639926897]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dT(self):
        test = self.phase.dA_dT(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[-5.08125375441449, -6.661433473344175]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dAj_dCik(self):
        test = self.phase.dAj_dCik(1700.0, 100000.0, [[1.0], [1.0], [1.0]], 0, 1, 0)
        answer = 0.0
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dC(self):
        test = self.phase.dGamma_i_dC(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[[[-0.0], [-0.0], [-0.0]], [[0.0], [0.0], [0.0]], [[0.0], [0.0], [0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dP(self):
        test = self.phase.dGamma_i_dP(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[-1.6588812253097058e-06, 9.157611533007939e-07, 7.43120072008912e-07]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dPhi(self):
        test = self.phase.dGamma_i_dPhi(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[[-0.0, 0.19179386193742448, 0.0], [0.0, -0.19179386193742448, 0.19553846400640473], [0.0, 0.0, -0.19553846400640473]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dT(self):
        test = self.phase.dGamma_i_dT(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[4.315410192723747e-05, 1.1505250124241901e-05, -5.465935205147937e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dC(self):
        test = self.phase.dGamma_ik_dC(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33], 1)
        answer =[[0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dP(self):
        test = self.phase.dGamma_ik_dP(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[[-1.6588812253097058e-06], [9.157611533007939e-07], [7.43120072008912e-07]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dPhi(self):
        test = self.phase.dGamma_ik_dPhi(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33], 1)
        answer =[[0.0, -0.19179386193742448, 0.19553846400640473]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dT(self):
        test = self.phase.dGamma_ik_dT(1700.0, 100000.0, [[1.0], [1.0], [1.0]], [0.33, 0.33, 0.33])
        answer =[[4.315410192723747e-05], [1.1505250124241901e-05], [-5.465935205147937e-05]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dMu_dC(self):
        test = self.phase.dMu_dC(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[[[0.0]], [[0.0]], [[0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dC(self):
        test = self.phase.ds_dC(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[[0.0], [0.0], [0.0]]
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
        answer =[[[-1.0], [1.0], [0.0]], [[0.0], [-1.0], [1.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_nu_m(self):
        test = self.phase.nu_m()
        answer =[[[-1.0], [1.0], [0.0]], [[0.0], [-1.0], [1.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_report(self):
        test = self.phase.report()
        answer = None
        assert(test == answer)

    def test_rho(self):
        test = self.phase.rho(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[33.450533005483145, 35.38916576209499, 36.33242393609047]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 100000.0, [[1.0], [1.0], [1.0]])
        answer =[2.5407332824429445, 2.5046174121236855, 2.4572701472376024]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_set_parameter(self):
        test = self.phase.set_parameter('R', 8.31442)
        answer = None
        assert(test == answer)

    def test_zero_C(self):
        test = self.phase.zero_C()
        answer =[[0.0], [0.0], [0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

