from __future__ import print_function, division
import numpy as np

from parent_class import EOS

class SoaveRedlichKwong(EOS):
    '''https://www.e-education.psu.edu/png520/m10_p5.html
    '''

    def __init__(self, Tc, pc, M, omega, fluid):

        super(SoaveRedlichKwong, self).__init__(M, fluid)
        self.acentric = omega
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a = 0.42748 * self.R**2 * self.T_crit**2 / self.p_crit
        self.b = 0.08664 * self.R * self.T_crit / self.p_crit

    def get_Z(self, T, p):

        kappa = 0.48508 + 1.55171*self.acentric - 0.15613*self.acentric**2
        alpha = (1 + kappa*(1-T/self.T_crit))**2

        Tr = T/self.T_crit
        alpha = (1 + kappa*(1 - Tr**0.5))**2

        A = alpha * self.a * p / self.R**2 / T**2
        B = self.b * p / self.R / T

        # Solve the cubic equation for compressibility factor z
        coeffs = [1, -1, A-B-B**2, -A*B]
        roots = np.roots(coeffs)

        real_roots = roots[np.isreal(roots)].real
        valid_roots = real_roots[real_roots > p*self.b/self.R/T]

        return valid_roots

