from __future__ import division, print_function
import numpy as np

class EOS(object):

    R = 8.3144598 # J/K/mol

    def __init__(self, Tc, pc, M):
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

    def __init__(self):
        pass

    def get_z(self, T, p):
        return 1.

    def get_rho(self, T, p):
        z = self.get_z(T, p)
        return p/z/self.R_sp/T

class RedlichKwong(EOS):
    '''https://www.e-education.psu.edu/png520/m10_p4.html
    '''

    def __init__(self, Tc, pc, M):

        #super(self.__class__, self).__init__(Tc, pc, M)
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a = 0.42748 * self.R**2 * self.T_crit**2.5 / self.p_crit
        self.b = 0.08664 * self.R * self.T_crit / self.p_crit

    def get_z(self, T, p):
        A = self.a * p / self.R**2 / T**2.5
        B = self.b * p / self.R / T

        # Solve the cubic equation for compressibility factor z
        coeffs = [-A*B, A-B-B**2, -1, 1]
        roots = np.roots(coeffs)

        real_roots = roots[np.isreal(roots)].real
        valid_roots = real_roots[real_roots > p*self.b/self.R/T]

        return valid_roots


class PengRobinson(EOS):
    '''https://www.e-education.psu.edu/png520/m11_p2.html
    '''

    def __init__(self, Tc, pc, M, omega):

        #super(self.__class__, self).__init__(Tc, pc, M)
        self.acentric = omega
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a = 0.45724 * self.R**2 * self.T_crit**2 / self.p_crit
        self.b = 0.07780 * self.R * self.T_crit / self.p_crit

    def get_z(self, T, p):

        kappa = 0.37464 + 1.54226*self.acentric - 0.26992*self.acentric**2
        alpha = (1 + kappa*(1-T/self.T_crit))**2

        A = self.a * p / self.R**2 / T**2.5
        B = self.b * p / self.R / T

        # Solve the cubic equation for compressibility factor z
        coeffs = [-(A*B-B**2-B**3), A-2*B-3*B**2, -(1-B), 1]
        roots = np.roots(coeffs)

        real_roots = roots[np.isreal(roots)].real
        valid_roots = real_roots[real_roots > p*self.b/self.R/T]

        return valid_roots
