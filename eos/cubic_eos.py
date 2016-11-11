from __future__ import print_function, division
import numpy as np

from cubic_parent import CubicEOS

class RK1(CubicEOS):
    '''https://www.e-education.psu.edu/png520/m11_p2.html
    '''

    def __init__(self, Tc, pc, M, fluid):

        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a0 = 0.42748 * self.R**2 * self.T_crit**2 / self.p_crit
        b1 = 0.08664 * self.R * self.T_crit / self.p_crit

        super(RK1, self).__init__(b1, 0., b1, 0., M, fluid)

    def get_a(self, T):
        Tr = T/self.T_crit
        return self.a0/Tr**0.5

    def get_diff_a_T(self, T):
        Tr = T/self.T_crit
        return -0.5*self.a0/T/Tr**0.5

    def get_double_diff_a_T(self, T):
        Tr = T/self.T_crit
        return 0.75*self.a0/T**2/Tr**0.5


class SRK1(CubicEOS):
    '''https://www.e-education.psu.edu/png520/m11_p2.html
    '''

    def __init__(self, Tc, pc, M, omega, fluid):

        self.acentric = omega
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a0 = 0.42748 * self.R**2 * self.T_crit**2 / self.p_crit
        b1 = 0.08664 * self.R * self.T_crit / self.p_crit
        self.kappa = 0.48508 + 1.55171*self.acentric - 0.15613*self.acentric**2

        super(SRK1, self).__init__(b1, 0., b1, 0., M, fluid)

    def get_a(self, T):
        Tr = T/self.T_crit
        alpha = (1 + self.kappa*(1 - Tr**0.5))**2
        return alpha * self.a0

    def get_diff_a_T(self, T):
        Tr = T/self.T_crit
        alpha0 = (1 + self.kappa*(1 - Tr**0.5))
        return -(self.a0*self.kappa/T)*Tr**0.5 * alpha0

    def get_double_diff_a_T(self, T):
        Tr = T/self.T_crit
        alpha0 = (1 + self.kappa*(1 - Tr**0.5))
        return (0.5*self.a0*self.kappa**2/T**2)*Tr + (0.5*self.a0*self.kappa/T**2)*Tr**0.5 * alpha0


class PR1(CubicEOS):
    '''https://www.e-education.psu.edu/png520/m11_p2.html
    '''

    def __init__(self, Tc, pc, M, omega, fluid):

        self.acentric = omega
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a0 = 0.45724 * self.R**2 * self.T_crit**2 / self.p_crit
        b1 = 0.07780 * self.R * self.T_crit / self.p_crit
        self.kappa = 0.37464 + 1.54226*self.acentric - 0.26992*self.acentric**2

        super(PR1, self).__init__(b1, 0., 2*b1, -b1**2, M, fluid)

    def get_a(self, T):
        Tr = T/self.T_crit
        alpha = (1 + self.kappa*(1 - Tr**0.5))**2
        return alpha * self.a0

    def get_diff_a_T(self, T):
        Tr = T/self.T_crit
        alpha0 = (1 + self.kappa*(1 - Tr**0.5))
        return -(self.a0*self.kappa/T)*Tr**0.5 * alpha0

    def get_double_diff_a_T(self, T):
        Tr = T/self.T_crit
        alpha0 = (1 + self.kappa*(1 - Tr**0.5))
        return (0.5*self.a0*self.kappa**2/T**2)*Tr + (0.5*self.a0*self.kappa/T**2)*Tr**0.5 * alpha0
