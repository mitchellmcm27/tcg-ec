import numpy as np
from numpy import array
from thermocodegen.testing import is_float_list, allclose_float_list
import py_fo_sio2_poly_linear_rxns as db
import pytest
class Testfo_sio2_poly_linear_rxns:
    phase = db.fo_sio2_poly_linear_rxns()

    def test_name(self):
        test = self.phase.name()
        answer = 'fo_sio2_poly_linear_rxns'
        assert(test == answer)

    def test_A(self):
        test = self.phase.A(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[-10117.936057890765, -3178.920714895008, 9904.611466794973]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_C_to_X(self):
        test = self.phase.C_to_X([[1.0], [1.0], [1.0], [0.25, 0.75]],)
        answer =[[1.0], [1.0], [1.0], [0.2807131061141594, 0.7192868938858406]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Cp(self):
        test = self.phase.Cp(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[1.3500876578848255, 1.3612319923012914, 1.2255008499285884, 1.7832122092897638]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_i(self):
        test = self.phase.Gamma_i(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[0.0018176021317915143, 0.0005710663751114755, -0.0, -0.0023886685069029897]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Gamma_ik(self):
        test = self.phase.Gamma_ik(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[[0.0018176021317915143], [0.0005710663751114755], [-0.0], [-0.00017089617559871675, -0.002217772331304273]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_M(self):
        test = self.phase.M()
        answer =[[140.69332], [100.388815], [60.08431], [120.16862, 140.69332]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_Mu(self):
        test = self.phase.Mu(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[[-2566345.4730803007], [-1823995.229304283], [-1075500.468689571], [-2170810.1603127318, -2556227.53702241]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_X_to_C(self):
        test = self.phase.X_to_C([[1.0], [1.0], [1.0], [0.2807131061141594, 0.7192868938858406]],)
        answer =[[1.0], [1.0], [1.0], [0.24999999999999994, 0.75]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_alpha(self):
        test = self.phase.alpha(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[5.1433326723180724e-05, 4.3568139995068164e-05, 3.178481588845701e-06, 7.696919774517023e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_beta(self):
        test = self.phase.beta(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[7.479182638846417e-07, 7.146031185909278e-07, 1.0960040772775672e-06, 3.0235423072759076e-06]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dC(self):
        test = self.phase.dA_dC(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[[[0.0], [0.0], [0.0], [49096.038333797514, -16365.346111265835]], [[0.0], [0.0], [0.0], [-6902.362153073569, 2300.7873843578573]], [[0.0], [0.0], [0.0], [-62900.76263994465, 20966.92087998155]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dP(self):
        test = self.phase.dA_dP(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[-0.12094378017311147, -0.37658754613824197, 0.15597746832738046]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dA_dT(self):
        test = self.phase.dA_dT(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[40.93604229062976, 21.943825383164395, 6.233927090809999]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dAj_dCik(self):
        test = self.phase.dAj_dCik(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], 0, 3, 0)
        answer = 49096.038333797514
        assert(np.isclose(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dC(self):
        test = self.phase.dGamma_i_dC(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[[[-0.0], [-0.0], [-0.0], [-0.00881969044155346, 0.0029398968138511527]], [[-0.0], [-0.0], [-0.0], [0.0012399513193245989, -0.0004133171064415332]], [[-0.0], [-0.0], [-0.0], [0.0, -0.0]], [[0.0], [0.0], [0.0], [0.007579739122228861, -0.0025265797074096197]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dP(self):
        test = self.phase.dGamma_i_dP(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[2.172653310040764e-08, 6.765078596569974e-08, -0.0, -8.937731906610737e-08]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dPhi(self):
        test = self.phase.dGamma_i_dPhi(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[[-0.0, 0.0, 0.0, 0.01817602131791514], [0.0, -0.0, 0.0, 0.005710663751114754], [0.0, 0.0, -0.01779280162832551, -0.0], [0.0, 0.0, 0.01779280162832551, -0.023886685069029894]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_i_dT(self):
        test = self.phase.dGamma_i_dT(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[-7.322369365563642e-06, -3.9321436649182625e-06, -0.0, 1.1254513030481905e-05]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dC(self):
        test = self.phase.dGamma_ik_dC(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1], 3)
        answer =[[-0.0003710653395759687, 0.00012368844652532295], [0.00795080446180483, -0.0026502681539349424]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dP(self):
        test = self.phase.dGamma_ik_dP(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[[2.172653310040764e-08], [6.765078596569974e-08], [-0.0], [-2.024503823312763e-08, -6.913228083297975e-08]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dPhi(self):
        test = self.phase.dGamma_ik_dPhi(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1], 3)
        answer =[[0.0, 0.0, 0.01779280162832551, -0.0017089617559871673], [0.0, 0.0, 0.0, -0.022177723313042727]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dGamma_ik_dT(self):
        test = self.phase.dGamma_ik_dT(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]], [0.6, 0.3, 0.0, 0.1])
        answer =[[-7.322369365563642e-06], [-3.9321436649182625e-06], [-0.0], [1.176725409735562e-06, 1.0077787620746343e-05]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_dMu_dC(self):
        test = self.phase.dMu_dC(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[[[0.0]], [[0.0]], [[0.0]], [[125801.5252798893, -41933.8417599631], [-49096.038333797514, 16365.346111265835]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_ds_dC(self):
        test = self.phase.ds_dC(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[[0.0], [0.0], [0.0], [2.7049759741363446, 2.94223079600028]]
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
        answer =[[[-1.0], [0.0], [0.0], [0.0, 1.0]], [[0.0], [-1.0], [0.0], [0.25, 0.5]], [[0.0], [0.0], [-1.0], [0.5, 0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_nu_m(self):
        test = self.phase.nu_m()
        answer =[[[-1.0], [0.0], [0.0], [0.0, 1.0]], [[0.0], [-1.0], [0.0], [0.2992579900460027, 0.7007420099539974]], [[0.0], [0.0], [-1.0], [1.0, 0.0]]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_report(self):
        test = self.phase.report()
        answer = None
        assert(test == answer)

    def test_rho(self):
        test = self.phase.rho(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[30.458046154746736, 30.55861099310908, 21.93493256414871, 27.764583855156445]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_s(self):
        test = self.phase.s(1700.0, 1000.0, [[1.0], [1.0], [1.0], [0.25, 0.75]])
        answer =[2.6512756103980992, 2.652644998492179, 2.6012229937579683, 2.882920591472729]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

    def test_set_parameter(self):
        test = self.phase.set_parameter('R', 8.314)
        answer = None
        assert(test == answer)

    def test_zero_C(self):
        test = self.phase.zero_C()
        answer =[[0.0], [0.0], [0.0], [0.0, 0.0]]
        assert(allclose_float_list(test,answer,rtol=1e-05,atol=1e-08))

