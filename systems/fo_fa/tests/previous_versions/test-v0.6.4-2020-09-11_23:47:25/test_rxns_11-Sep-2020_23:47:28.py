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
        answer =[-18328.996493542567, -8162.250289141666]
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
        answer =[0.8062166068802339, -0.8062166068802339]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_ik(self):
        test = self.phase.Gamma_ik(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[[0.5578122268749744, 0.2484043800052595], [-0.5578122268749744, -0.2484043800052595]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_M(self):
        test = self.phase.M()
        answer =[[140.69332, 203.77312], [140.69332, 203.77312]]
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
        answer =[[[1118.9614607218273, -10070.653146496446], [-7425.412798522721, 11138.11919778408]], [[-14586.134595789852, 131275.21136210865], [16132.231286244789, -24198.34692936718]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dP(self):
        test = self.phase.dA_dP(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[-0.3613726772634793, -0.5374277759668677]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dT(self):
        test = self.phase.dA_dT(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[40.73002905867378, 52.6480504673267]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dAj_dCik(self):
        test = self.phase.dAj_dCik(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], 0, 1, 0)
        answer = -7425.412798522721
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dC(self):
        test = self.phase.dGamma_i_dC(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[[[0.4098507978235312, -3.68865718041178], [-0.26497739859045305, 0.3974660978856795]], [[-0.4098507978235312, 3.68865718041178], [0.26497739859045305, -0.3974660978856795]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dP(self):
        test = self.phase.dGamma_i_dP(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[2.7353482363821148e-05, -2.7353482363821148e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dPhi(self):
        test = self.phase.dGamma_i_dPhi(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[[-0.0, 0.671847172400195], [0.0, -0.671847172400195]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dT(self):
        test = self.phase.dGamma_i_dT(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[-0.0030370824116391594, 0.0030370824116391594]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dC(self):
        test = self.phase.dGamma_ik_dC(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 1)
        answer =[[-0.2259799684106626, 0.3389699526159939], [0.49095736700111564, -0.7364360505016734]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dP(self):
        test = self.phase.dGamma_ik_dP(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[[1.0997770549365876e-05, 1.635571181445527e-05], [-1.0997770549365876e-05, -1.635571181445527e-05]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dPhi(self):
        test = self.phase.dGamma_ik_dPhi(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6], 1)
        answer =[[0.0, -0.46484352239581206], [0.0, -0.20700365000438292]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dT(self):
        test = self.phase.dGamma_ik_dT(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]], [0.4, 0.6])
        answer =[[-0.0013746602610431975, -0.0016624221505959621], [0.0013746602610431975, 0.0016624221505959621]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dMu_dC(self):
        test = self.phase.dMu_dC(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[[[1118.9614607218273, -10070.653146496446], [-14586.134595789852, 131275.21136210865]], [[7425.412798522721, -11138.11919778408], [-16132.231286244789, 24198.34692936718]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dC(self):
        test = self.phase.ds_dC(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
        answer =[[2.655641037934999, 2.329758657117597], [2.9451382049881367, 2.588131437588576]]
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
        answer =[[[-1.0, 0.0], [1.0, 0.0]], [[0.0, -1.0], [0.0, 1.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_nu_m(self):
        test = self.phase.nu_m()
        answer =[[[-1.0, 0.0], [1.0, 0.0]], [[0.0, -1.0], [0.0, 1.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_report(self):
        test = self.phase.report()
        answer = None
        assert(test == answer)

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, [[0.9, 0.1], [0.6, 0.4]])
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

