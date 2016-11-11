from __future__ import print_function, division
import numpy as np

from parent_class import EOS

class PengRobinson(EOS):
    '''https://www.e-education.psu.edu/png520/m11_p2.html
    '''

    def __init__(self, Tc, pc, M, omega, fluid):

        super(PengRobinson, self).__init__(M, fluid)
        self.acentric = omega
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

        self.a = 0.45724 * self.R**2 * self.T_crit**2 / self.p_crit
        self.b = 0.07780 * self.R * self.T_crit / self.p_crit

    def get_Z(self, T, p):

        kappa = 0.37464 + 1.54226*self.acentric - 0.26992*self.acentric**2

        Tr = T/self.T_crit
        alpha = (1 + kappa*(1 - Tr**0.5))**2

        A = alpha * self.a * p / self.R**2 / T**2
        B = self.b * p / self.R / T

        # Solve the cubic equation for compressibility factor z
        coeffs = [1, -(1-B), A-2*B-3*B**2, -(A*B-B**2-B**3)]
        #print(coeffs)
        roots = np.roots(coeffs)

        real_roots = roots[np.isreal(roots)].real
        valid_roots = real_roots[real_roots > p*self.b/self.R/T]

        return valid_roots

