from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

from customize_matplotlib import *
customize_matplotlib()

T = 40.0
show_legend = True
# show_legend = False
# T = 60.0
# T = 80.0
# T = 100.0


rho = np.loadtxt('data/rho_{}.txt'.format(T))
c = np.loadtxt('data/c_{}.txt'.format(T))

labels = ['Ideal gas', 'RK', 'SRK', 'PR', 'LK', 'SW']
#styles = [':ko', ':ks', ':k*', ':k<', ':kd', '-k']
styles = ['-.k', '--k', ':xk', ':+k', ':*k', '-k']

fig_r, ax_r = plt.subplots()
fig_c, ax_c = plt.subplots()

ax_r.set_xlabel('Pressure [Mpa]')
ax_r.set_ylabel(r'Density [kg/m$^3$]')
ax_c.set_xlabel('Pressure [Mpa]')
ax_c.set_ylabel('Speed of sound [m/s]')

for i, label in enumerate(labels):
    ax_r.plot(rho[:,0]/1e6, rho[:,i+1], styles[i], label=label, fillstyle='none')
    ax_c.plot(c[:,0]/1e6, c[:,i+1], styles[i], label=label, fillstyle='none')

if show_legend:
    ax_r.legend(loc='best')
    ax_c.legend(loc='best')

fig_r.savefig('plots/rho_{}.pdf'.format(T))
fig_c.savefig('plots/c_{}.pdf'.format(T))

plt.show()

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
