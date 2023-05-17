import io
import numpy as np
from contextlib import redirect_stdout, redirect_stderr
import time

from ..base import *
from .base import BasePDReactiveODE
from scipy.integrate import solve_ivp, BDF, OdeSolution
from scipy.optimize import OptimizeResult as OdeResult

SCIPY_AVAIL = True
MESSAGES = {
    0: "The solver successfully reached the end of the integration interval.",
    1: "A termination event occurred.",
    2: "The solver failed to converge."
}

class ScipyPDReactiveODE(BasePDReactiveODE):

    def solve(self,T,p,mi0,Cik0,end,**kwargs):

        max_steps = kwargs.pop("max_steps", np.inf)

        super().set_initial_params(T,p,mi0,Cik0,**kwargs)

        method = kwargs.get('method', 'BDF')
        rtol   = kwargs.get('rtol', 1.e-5)
        atol   = kwargs.get('atol', 1.e-9)

        u0 = np.concatenate((mi0,Cik0))

        #try:
        #so = io.StringIO()
        #se = io.StringIO()
        #with redirect_stdout(so), redirect_stderr(se):
        tic = time.perf_counter()
        if method == 'BDF_mcm':
            t0 = 0.
            solver = BDF(self.rhs, t0, u0, float(end), rtol=rtol, atol=atol, jac=self.jac)
            ts = [t0]
            ys = [u0]
            interpolants = []
            status = None
            num_steps = 0
            while status is None:
                message = solver.step()
          
                if solver.status == 'finished':
                    status = 0
                elif solver.status == 'failed':
                    status = -1
                    break

                if num_steps == max_steps:
                    status = 2
                
                t_old = solver.t_old
                t = solver.t
                y = solver.y
                output = solver.dense_output()
                interpolants.append(output)
                ts.append(t)
                ys.append(y)
                num_steps = num_steps + 1
                
            message = MESSAGES.get(status, message)
            ts = np.array(ts)
            ys = np.vstack(ys).T
            ode_solution = OdeSolution(ts, interpolants)
            sol = OdeResult(t=ts, y=ys, sol=ode_solution,
                    nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
                    status=status, message=message, success=status >= 0)
        else:
          sol = solve_ivp(self.rhs, [0,end], u0, dense_output=True, method=method, rtol=rtol, atol=atol, jac=self.jac)
        
        self.sol = sol
        toc = time.perf_counter()
        self.stime = toc-tic
        #self.stdout = so.getvalue()
        #self.stderr = se.getvalue()
        #flag = self.rxn.check_coder_error()
          #if flag!=0:
          #  self.sol = None
          #  self.excstr = repr(flag)
          #  self.rxn.reset_coder_error()
        # except Exception as e:
        #   self.sol = None
        #  self.excstr = repr(e)

