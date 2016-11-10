from __future__ import print_function, division
import numpy as np
from math import log

from parent_class import EOS

class RedlichKwong(EOS):
    '''https://www.e-education.psu.edu/png520/m10_p4.html
    '''

    def __init__(self, Tc, pc, M, fluid):

        super(RedlichKwong, self).__init__(M, fluid)
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a = 0.42748 * self.R**2 * self.T_crit**2.5 / self.p_crit
        self.b = 0.08664 * self.R * self.T_crit / self.p_crit

    # Coefficient A
    def get_A(self, T, p):
        return self.a * p / self.R**2 / T**2.5

    # dA/dT at const. p
    def get_pdiff_A_T_p(self, T, p):
        return -2.5 * self.a * p / self.R**2 / T**3.5

    # dA/dp at const. T
    def get_pdiff_A_p_T(self, T, p):
        return self.a / self.R**2 / T**2.5

    # dA/dT at const. v
    #def get_pdiff_A_T_v(self, T, p):
    #    return -2.5 * self.a * p / self.R**2 / T**3.5

    # Coefficient B
    def get_B(self, T, p):
        return self.b * p / self.R / T

    # dB/dT at const. p
    def get_pdiff_B_T_p(self, T, p):
        return -1 * self.b * p / self.R / T**2

    # dB/dp at const. T
    def get_pdiff_B_p_T(self, T, p):
        return self.b / self.R / T

    def get_Z(self, T, p):
        #A = self.a * p / self.R**2 / T**2.5
        #B = self.b * p / self.R / T
        A = self.get_A(T, p)
        B = self.get_B(T, p)

        # Solve the cubic equation for compressibility factor z
        # Z^3 - Z^2 + (A-B-B**2)*Z - A*B = 0
        coeffs = [1, -1, A-B-B**2, -A*B]
        roots = np.roots(coeffs)

        real_roots = roots[np.isreal(roots)].real
        valid_roots = real_roots[real_roots > p*self.b/self.R/T]

        return valid_roots

    # dZ/dT at const. p
    def get_pdiff_Z_T_p(self, T, p, **kwargs):
        A = self.get_A(T, p)
        B = self.get_B(T, p)
        A1 = self.get_pdiff_A_T_p(T, p)
        B1 = self.get_pdiff_B_T_p(T, p)

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = get_Z(self, T, p)

        return (A*B1 + A1*B - (A1-B1-2*B*B1)) * Z / (3*Z**2 - 2*Z + (A-B-B**2))

    # dZ/dp at const. T
    def get_pdiff_Z_p_T(self, T, p, **kwargs):
        A = self.get_A(T, p)
        B = self.get_B(T, p)
        A1 = self.get_pdiff_A_p_T(T, p)
        B1 = self.get_pdiff_B_p_T(T, p)

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = get_Z(self, T, p)

        return (A*B1 + A1*B - (A1-B1-2*B*B1)) * Z / (3*Z**2 - 2*Z + (A-B-B**2))

    def get_departure_cp(self, T, p, **kwargs):

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        Z1 = self.get_pdiff_Z_T_p(T, p, Z=Z)

        h = self.b*p/Z/self.R/T

        return self.R*(Z + T*Z1 -1) - 3*self.a/2/self.b * (-log(1+h)/2/T**1.5 + 1/T**0.5/(1+h) * self.b*p/self.R * (-Z-T*Z1)/Z**2/T**2)


    def get_departure_cv(self, T, p, **kwargs):

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        h = self.b*p/Z/self.R/T

        return -3*self.a/2/self.b * (-log(1+h)/2/T**1.5)

    def get_isothermal_compressibility(self, T, p, **kwargs):

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = get_Z(self, T, p)

        Z1 = self.get_pdiff_Z_p_T(T, p, Z=Z)

        return 1/p - 1/Z * Z1

