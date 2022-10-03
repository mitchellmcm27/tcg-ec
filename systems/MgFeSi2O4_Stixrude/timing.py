import numpy as np
import timeit
import py_Mg2SiO4_stixrude as pms

rxn = pms.Mg2SiO4_stixrude()
rxn.report()

T = 2000.
P = 250000.
C = np.ones((3,1))
N=10000

total_time = 0. 
for method in ['A', 'rho', 's', 'alpha', 'Cp', 'beta']:
    time = timeit.timeit('rxn.{}(T,P,C)'.format(method),globals=globals(),number=N)/N
    print('{}:\t {}'.format(method,time))
    total_time += time

print('\nTotal time = {}'.format(total_time))




