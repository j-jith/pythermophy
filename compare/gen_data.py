from __future__ import print_function, division
import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

import sys
sys.path.append('../')
import eos

### Mean temperature of fluid
T_0 = np.array([40, 60, 80, 100]) + 273.15 # [K]
### Mean pressure of fluid
p_0 = np.linspace(0.1e6, 20e6, 50) # [Pa]

# CO2 properties
M = 44.01e-3 # kg/mol
pc = 7.38e6 # Pa
Tc = 31.1 + 273.15 # K
omega = 0.228 # acentric factor

ig = eos.IdealGas(M, 'CO2')
rk = eos.RK1(Tc, pc, M, 'CO2')
srk = eos.SRK1(Tc, pc, M, omega, 'CO2')
pr = eos.PR1(Tc, pc, M, omega, 'CO2')
lk = eos.LK(Tc, pc, M, omega, 'CO2')

eqs = [ig, rk, srk, pr, lk]
eqs_labels = ['Ideal gas', 'RK', 'SRK', 'PR', 'LK']
eqs_colors = ['b', 'g', 'r', 'c', 'y']


rho = np.empty((len(p_0), len(eqs)+1))
c = np.empty((len(p_0), len(eqs)+1))

for T in T_0:
    for j, eq in enumerate(eqs):
        for i, p in enumerate(p_0):
            rho[i, j] = eq.get_rho(T, p)
            c[i, j] = eq.get_speed_of_sound(T, p)

    rho[:, j+1] = PropsSI('DMASS', 'T', T, 'P', p_0, 'CO2')
    c[:, j+1] = PropsSI('SPEED_OF_SOUND', 'T', T, 'P', p_0, 'CO2')

    header = 'T = {} C (+273.15 K)\np, Ideal gas, RK, SRK, PR, LK, SW'.format(T-273.15)
    np.savetxt('data/rho_{}.txt'.format(T-273.15), np.hstack([p_0[:, None], rho]), header=header)
    np.savetxt('data/c_{}.txt'.format(T-273.15), np.hstack([p_0[:, None], c]), header=header)


# markers = ['o', 's', '*', '<', 'd']
# 
# fig, ax = plt.subplots(nrows=2, ncols=2)
# 
# n = 0
# for i in range(2):
#     for j in range(2):
#         for k, eq in enumerate(eqs):
#             rho = []
#             for p in p_0:
#                 rho.append(eq.get_rho(T_0[n], p))
#             ax[i][j].plot(p_0/1e6, rho, ':k'+markers[k], fillstyle='none', label=eqs_labels[k])
# 
#         rho_sw = PropsSI('DMASS', 'T', T_0[n], 'P', p_0, 'CO2')
#         ax[i][j].plot(p_0/1e6, rho_sw, '-k', label='S & W')
# 
#         ax[i][j].set_xlabel('Pressure [MPa]')
#         ax[i][j].set_ylabel(r'Density [kg/m$^3$]')
#         #ax[i][j].set_title(r'{} $^{{\circ}}$C'.format(T_0[n]-273.15))
#         ax[i][j].text(.5,.9, r'{} $^{{\circ}}$C'.format(T_0[n]-273.15), horizontalalignment='center', transform=ax[i][j].transAxes)
#         n+=1
# 
# plt.show()
# 
# 
# # for i, eq in enumerate(eqs):
# #     for j, T in enumerate(T_0):
# #     #T = T_0[0]
# #         prop = []
# #         for p in p_0:
# #             #prop.append(eq.get_rho(T, p))
# #             prop.append(eq.get_speed_of_sound(T, p))
# # 
# #         ax.plot(p_0, prop, '-'+eqs_colors[i]+T_styles[j], label=eqs_labels[i])
# # 
# # for j, T in enumerate(T_0):
# #     #sw = PropsSI('DMASS', 'T', T, 'P', p_0, 'CO2')
# #     sw = PropsSI('SPEED_OF_SOUND', 'T', T, 'P', p_0, 'CO2')
# #     ax.plot(p_0, sw, '-k'+T_styles[j], label='S & W')
# # 
# # ax.legend(loc='best')
# # plt.show()
