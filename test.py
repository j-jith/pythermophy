from __future__ import print_function, division
import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

import eos

### Mean temperature of fluid
T_0 = 328.79 # [K]
### Mean pressure of fluid
p_0 = np.linspace(1e6, 20e6, 10) # [Pa]

# CO2 properties
M = 44.01e-3 # kg/mol
pc = 7.38e6 # Pa
Tc = 31.1 + 273.15 # K
omega = 0.228 # acentric factor

ig = eos.IdealGas(M, 'CO2')
rk = eos.RK(Tc, pc, M, 'CO2')
srk = eos.SRK(Tc, pc, M, omega, 'CO2')
pr = eos.PR(Tc, pc, M, omega, 'CO2')
lk = eos.LK(Tc, pc, M, omega, 'CO2')

i = 9

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
#print('cp: ', rk.get_cp(T_0, p_0[i], Z=Z)/M)
#print('cv: ', rk.get_cv(T_0, p_0[i], Z=Z)/M)
#print('beta: ', rk.get_isothermal_compressibility(T_0, p_0[i], Z=Z))
#print('gamma: ', rk.get_adiabatic_index(T_0, p_0[i]))
#print('c: ', rk.get_speed_of_sound(T_0, p_0[i], Z=Z))

print('\nSpan-Wagner')
print('rho: ', PropsSI('DMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('cp: ', PropsSI('CPMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('cv: ', PropsSI('CVMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('beta: ', PropsSI('ISOTHERMAL_COMPRESSIBILITY', 'T', T_0, 'P', p_0[i], 'CO2'))
print('gamma: ', PropsSI('CPMASS', 'T', T_0, 'P', p_0[i], 'CO2')/PropsSI('CVMASS', 'T', T_0, 'P', p_0[i], 'CO2'))
print('c: ', PropsSI('SPEED_OF_SOUND', 'T', T_0, 'P', p_0[i], 'CO2'))

# rho_0 = PropsSI('DMASS', 'T', T_0, 'P', p_0, 'CO2')
# 
# #print(rk.get_rho(T_0, p_0[-1]))
# #print(srk.get_z(T_0, p_0[-1]))
# #print(lk.get_rho(T_0, p_0[-1]))
# 
# #lk.get_z(T_0, p_0[-1])
# 
# 
# rho = []
# for i, p in enumerate(p_0):
#     rho.append([ig.get_rho(T_0, p), rk.get_rho(T_0, p)[0], srk.get_rho(T_0, p)[0], max(pr.get_rho(T_0, p)), lk.get_rho(T_0, p)[0], rho_0[i]])
# #print(rho)
# 
# rho = np.array(rho)
# 
# fig, ax = plt.subplots()
# 
# labels = ['Ideal gas', 'RK', 'SRK', 'PR', 'LK', 'SW']
# 
# for i in range(len(rho[0])):
#     ax.plot(p_0/1e6, rho[:, i], label=labels[i])
# 
# ax.set_xlabel('Pressure [MPa]')
# ax.set_ylabel(r'Density [kg/m$^3$]')
# ax.legend(loc='best')
# 
# plt.show()
