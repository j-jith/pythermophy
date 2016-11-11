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

pr = eos.PR(Tc, pc, M, omega, 'CO2')
pr1 = eos.PR1(Tc, pc, M, omega, 'CO2')


rk = eos.RK(Tc, pc, M, 'CO2')
rk1 = eos.RK1(Tc, pc, M, 'CO2')

srk = eos.SRK(Tc, pc, M, omega, 'CO2')
srk1 = eos.SRK1(Tc, pc, M, omega, 'CO2')

i = 5

print(pr.get_Z(T_0, p_0[i]))
print(pr1.get_Z(T_0, p_0[i]))

#print(rk.get_speed_of_sound(T_0, p_0[i]))
#print(rk1.get_speed_of_sound(T_0, p_0[i]))
#print(PropsSI('SPEED_OF_SOUND', 'T', T_0, 'P', p_0[i], 'CO2'))
