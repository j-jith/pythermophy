from __future__ import print_function, division
import numpy as np

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

rk = eos.RedlichKwong(Tc, pc, M)
pr = eos.PengRobinson(Tc, pc, M, omega)

print(rk.get_rho(T_0, p_0[-1]))
print(pr.get_rho(T_0, p_0[-1]))

#z_rk = np.empty(p_0.shape)
#for i in range(len(p_0)):
#    z_rk[i] = rk.get_z(T_0, p_0[i])

