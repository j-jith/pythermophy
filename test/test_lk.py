from __future__ import print_function, division
import numpy as np
import pythermophy as eos

### Mean temperature of fluid
#T_0 = 328.79 # [K]
T_0 = 100 + 273.15 # [K]
### Mean pressure of fluid
p_0 = np.linspace(1e6, 20e6, 10) # [Pa]

# CO2 properties
#M = 44.01e-3 # kg/mol
#pc = 7.38e6 # Pa
#Tc = 31.1 + 273.15 # K
#omega = 0.228 # acentric factor

fluid = eos.Fluid.init_from_file('fluids/CO2')
M = fluid.molar_mass
lk = eos.LK(fluid)

i = -1

# print('Cp')
# print(lk.get_cp(T_0, p_0[i], step=1e-2)/M)
# print(lk.get_cp(T_0, p_0[i], step=1e-3)/M)
# print(lk.get_cp(T_0, p_0[i], step=1e-4)/M)
# 
# print('Cv')
# print(lk.get_cv(T_0, p_0[i], step=1e-2)/M)
# print(lk.get_cv(T_0, p_0[i], step=1e-3)/M)
# print(lk.get_cv(T_0, p_0[i], step=1e-4)/M)
# 
# print('dZ/dp at const. T')
# print(lk.get_pdiff_Z_p_T(T_0, p_0[i], step=1e-2)/M)
# print(lk.get_pdiff_Z_p_T(T_0, p_0[i], step=1e-3)/M)
# print(lk.get_pdiff_Z_p_T(T_0, p_0[i], step=1e-4)/M)
