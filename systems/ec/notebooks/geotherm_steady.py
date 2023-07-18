# %%

import numpy as np
import matplotlib.pyplot as plt
from fd1d_heat_implicit import *

b1 = 10 + 273.15
b2 = 1330 + 273.15

def geotherm_steady(x,L,shortening):
    # x is 0 to 1
    k = 3.35
    A = 1.9e-6 # rhoH0
    b = 10.e3 * shortening # hr

    a1 = b1 + b**2/k*A
    a0 = b2 - a1 + b**2/k*A*np.exp(-L/b)
    Te = a0*x + a1 - b**2/k*A*np.exp(-L/b*x)

    dTdx = a0 + b*A*L/k*np.exp(-L/b*x)
    return Te, dTdx/L
# %%
