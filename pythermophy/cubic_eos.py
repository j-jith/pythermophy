from __future__ import print_function, division
import numpy as np

from cubic_parent import CubicEOS

class RedlichKwong(CubicEOS):

    def __init__(self, fluid):
        """
        Redlich-Kwong equation of state
        For details see: https://www.e-education.psu.edu/png520/m10_p4.html

        Parameters
        ----------
        fluid - A Fluid class object

        Returns
        -------
        An equation of state object
        """

        self.p_crit = fluid.p_crit # Pa
        self.T_crit = fluid.T_crit # K

        self.a0 = 0.42748 * self.R**2 * self.T_crit**2 / self.p_crit
        b1 = 0.08664 * self.R * self.T_crit / self.p_crit

        super(RedlichKwong, self).__init__(b1, 0., b1, 0., fluid)

    def get_a(self, T):
        Tr = T/self.T_crit
        return self.a0/Tr**0.5

    def get_diff_a_T(self, T):
        Tr = T/self.T_crit
        return -0.5*self.a0/T/Tr**0.5

    def get_double_diff_a_T(self, T):
        Tr = T/self.T_crit
        return 0.75*self.a0/T**2/Tr**0.5


class SoaveRedlichKwong(CubicEOS):

    def __init__(self, fluid):
        """
        Soave-Redlich-Kwong equation of state
        For details see: https://www.e-education.psu.edu/png520/m10_p5.html

        Parameters
        ----------
        fluid - A Fluid class object

        Returns
        -------
        An equation of state object
        """

        self.acentric = fluid.acentric
        self.p_crit = fluid.p_crit # Pa
        self.T_crit = fluid.T_crit # K

        self.a0 = 0.42748 * self.R**2 * self.T_crit**2 / self.p_crit
        b1 = 0.08664 * self.R * self.T_crit / self.p_crit
        self.kappa = 0.48508 + 1.55171*self.acentric - 0.15613*self.acentric**2

        super(SoaveRedlichKwong, self).__init__(b1, 0., b1, 0., fluid)

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


class PengRobinson(CubicEOS):

    def __init__(self, fluid):
        """
        Peng-Robinson equation of state
        For details see: https://www.e-education.psu.edu/png520/m11_p2.html

        Parameters
        ----------
        fluid - A Fluid class object

        Returns
        -------
        An equation of state object
        """

        self.acentric = fluid.acentric
        self.p_crit = fluid.p_crit # Pa
        self.T_crit = fluid.T_crit # K

        self.a0 = 0.45724 * self.R**2 * self.T_crit**2 / self.p_crit
        b1 = 0.07780 * self.R * self.T_crit / self.p_crit
        self.kappa = 0.37464 + 1.54226*self.acentric - 0.26992*self.acentric**2

        super(PengRobinson, self).__init__(b1, 0., 2*b1, -b1**2, fluid)

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
