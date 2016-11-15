from __future__ import print_function, division
import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

import sys
sys.path.append('../')
import pythermophy as eos

### Mean temperature of fluid
T_0 = 328.79 # [K]
### Mean pressure of fluid
p_0 = np.linspace(1e6, 20e6, 10) # [Pa]

# CO2 properties
#M = 44.01e-3 # kg/mol
#pc = 7.38e6 # Pa
#Tc = 31.1 + 273.15 # K
#omega = 0.228 # acentric factor

fluid = eos.Fluid.init_from_file('../fluids/CO2')
M = fluid.molar_mass

ig = eos.IdealGas(fluid)
rk = eos.RK(fluid)
srk = eos.SRK(fluid)
pr = eos.PR(fluid)
lk = eos.LK(fluid)

i = -1

print('Ideal gas')
Z = ig.get_Z(T_0, p_0[i])
print('Z: ', Z)
print('rho: ', ig.get_rho(T_0, p_0[i]))
print('cp: ', ig.get_cp(T_0, p_0[i])/M)
print('cv: ', ig.get_cv(T_0, p_0[i])/M)
print('beta: ', ig.get_isothermal_compressibility(T_0, p_0[i]))
print('gamma: ', ig.get_adiabatic_index(T_0, p_0[i]))
print('c: ', ig.get_speed_of_sound(T_0, p_0[i]))

print('\nRK')
Z = rk.get_Z(T_0, p_0[i])
print('Z: ', Z)
print('rho: ', rk.get_rho(T_0, p_0[i], Z=Z))
print('cp: ', rk.get_cp(T_0, p_0[i], Z=Z)/M)
print('cv: ', rk.get_cv(T_0, p_0[i], Z=Z)/M)
print('beta: ', rk.get_isothermal_compressibility(T_0, p_0[i], Z=Z))
print('gamma: ', rk.get_adiabatic_index(T_0, p_0[i]))
print('c: ', rk.get_speed_of_sound(T_0, p_0[i], Z=Z))

print('\nSRK')
Z = srk.get_Z(T_0, p_0[i])
print('Z: ', Z)
print('rho: ', srk.get_rho(T_0, p_0[i], Z=Z))
print('cp: ', srk.get_cp(T_0, p_0[i], Z=Z)/M)
print('cv: ', srk.get_cv(T_0, p_0[i], Z=Z)/M)
print('beta: ', srk.get_isothermal_compressibility(T_0, p_0[i], Z=Z))
print('gamma: ', srk.get_adiabatic_index(T_0, p_0[i]))
print('c: ', srk.get_speed_of_sound(T_0, p_0[i], Z=Z))

print('\nPR')
Z = pr.get_Z(T_0, p_0[i])
print('Z: ', Z)
print('rho: ', pr.get_rho(T_0, p_0[i], Z=Z))
print('cp: ', pr.get_cp(T_0, p_0[i], Z=Z)/M)
print('cv: ', pr.get_cv(T_0, p_0[i], Z=Z)/M)
print('beta: ', pr.get_isothermal_compressibility(T_0, p_0[i], Z=Z))
print('gamma: ', pr.get_adiabatic_index(T_0, p_0[i]))
print('c: ', pr.get_speed_of_sound(T_0, p_0[i], Z=Z))

print('\nLK')
Z = lk.get_Z(T_0, p_0[i])
print('Z: ', Z)
print('rho: ', lk.get_rho(T_0, p_0[i], Z=Z))
# print('cp: ', lk.get_cp(T_0, p_0[i], Z=Z)/M)
# print('cv: ', lk.get_cv(T_0, p_0[i], Z=Z)/M)
# print('beta: ', lk.get_isothermal_compressibility(T_0, p_0[i], Z=Z))
# print('gamma: ', lk.get_adiabatic_index(T_0, p_0[i]))
# print('c: ', lk.get_speed_of_sound(T_0, p_0[i], Z=Z))
print('cp: ', lk.get_cp(T_0, p_0[i])/M)
print('cv: ', lk.get_cv(T_0, p_0[i])/M)
print('beta: ', lk.get_isothermal_compressibility(T_0, p_0[i]))
print('gamma: ', lk.get_adiabatic_index(T_0, p_0[i]))
print('c: ', lk.get_speed_of_sound(T_0, p_0[i]))

print('\nSpan-Wagner')
print('rho: ', PropsSI('DMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('cp: ', PropsSI('CPMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('cv: ', PropsSI('CVMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('beta: ', PropsSI('ISOTHERMAL_COMPRESSIBILITY', 'T', T_0, 'P', p_0[i], 'CO2'))
print('gamma: ', PropsSI('CPMASS', 'T', T_0, 'P', p_0[i], 'CO2')/PropsSI('CVMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('c: ', PropsSI('SPEED_OF_SOUND', 'T', T_0, 'P', p_0[i], 'CO2'))

