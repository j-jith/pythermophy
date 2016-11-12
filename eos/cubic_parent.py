from __future__ import print_function, division
import numpy as np
from math import log

from parent_class import EOS

class CubicEOS(EOS):

    def __init__(self, b, c, d, e, fluid):

        self.b = b; self.c = c
        self.d = d; self.e = e

        super(CubicEOS, self).__init__(fluid)

    def get_a(self, T):
        """
        Temperature dependent coefficient "a". This function
        has to be redefined in the child class
        """
        pass

    def get_diff_a_T(self, T):
        """
        Derivative of coefficient "a" wrt temperature. This function
        has to be redefined in the child class
        """
        pass

    def get_double_diff_a_T(self, T):
        """
        Second derivative of coefficient "a" wrt temperature. This function
        has to be redefined in the child class
        """
        pass

    # Coefficients of the cubic polynomial for Z: A, B, C
    # Z^3 + A*Z^2 + B*Z + C = 0
    def get_A(self, T, p):
        return (-self.R*T - self.b*p + self.c*p + self.d*p)/(self.R*T)

    def get_B(self, T, p):
        a = self.get_a(T)
        return (-self.R*T*self.d*p + a*p - self.b*self.d*p**2 + self.c*self.d*p**2 + self.e*p**2)/(self.R**2*T**2)

    def get_C(self, T, p):
        a = self.get_a(T)
        return (-self.R*T*self.e*p**2 - a*self.b*p**2 + a*self.c*p**2 - self.b*self.e*p**3 + self.c*self.e*p**3)/(self.R**3*T**3)

    # Partial derivatives of A, B and C wrt. T at constant p
    def get_pdiff_A_T_p(self, T, p):
        return p*(self.b - self.c - self.d)/(self.R*T**2)

    def get_pdiff_B_T_p(self, T, p):
        a = self.get_a(T)
        dadT = self.get_diff_a_T(T)
        return p*(self.R*T*self.d + T*dadT + 2*self.b*self.d*p - 2*self.c*self.d*p - 2*self.e*p - 2*a)/(self.R**2*T**3)

    def get_pdiff_C_T_p(self, T, p):
        a = self.get_a(T)
        dadT = self.get_diff_a_T(T)
        return p**2*(2*self.R*T*self.e + T*(-self.b + self.c)*dadT + 3*self.b*self.e*p + 3*self.b*a - 3*self.c*self.e*p - 3*self.c*a)/(self.R**3*T**4)

    # Partial derivatives of A, B and C wrt. p at constant T
    def get_pdiff_A_p_T(self, T, p):
        return (-self.b + self.c + self.d)/(self.R*T)

    def get_pdiff_B_p_T(self, T, p):
        a = self.get_a(T)
        return (-self.R*T*self.d - 2*self.b*self.d*p + 2*self.c*self.d*p + 2*self.e*p + a)/(self.R**2*T**2)

    def get_pdiff_C_p_T(self, T, p):
        a = self.get_a(T)
        return p*(-2*self.R*T*self.e - 3*self.b*self.e*p - 2*self.b*a + 3*self.c*self.e*p + 2*self.c*a)/(self.R**3*T**3)


    # Solve the cubic equation for compressiblilty factor Z
    def get_Z(self, T, p):
        """
        Get the compressibility factor for a real gas
        Parameters
        ----------
        T - Temperature [K]
        p - Pressure [Pa]

        Returns
        -------
        Compressibility factor [dimensionless]
        """

        # a = self.get_a(T)
        # coeffs = [
        #         1,
        #         (-self.R*T - self.b*p + self.c*p + self.d*p)/(self.R*T),
        #         (-self.R*T*self.d*p + a*p - self.b*self.d*p**2 + self.c*self.d*p**2 + self.e*p**2)/(self.R**2*T**2),
        #         (-self.R*T*self.e*p**2 - a*self.b*p**2 + a*self.c*p**2 - self.b*self.e*p**3 + self.c*self.e*p**3)/(self.R**3*T**3)
        # ]
        coeffs = [1, self.get_A(T, p), self.get_B(T, p), self.get_C(T, p)]
        #print(coeffs)

        roots = np.roots(coeffs)
        real_roots = roots[np.isreal(roots)].real

        if len(real_roots) == 1:
            real_roots = real_roots[0]

        return real_roots

    # Partial derivative of Z wrt T at constant p
    def get_pdiff_Z_T_p(self, T, p, **kwargs):
        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        A = self.get_A(T, p)
        B = self.get_B(T, p)
        C = self.get_C(T, p)

        A1 = self.get_pdiff_A_T_p(T, p)
        B1 = self.get_pdiff_B_T_p(T, p)
        C1 = self.get_pdiff_C_T_p(T, p)

        return -(Z**2*A1 + Z*B1 + C1)/(2*A*Z + B + 3*Z**2)

    # Partial derivative of Z wrt p at constant T
    def get_pdiff_Z_p_T(self, T, p, **kwargs):
        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        A = self.get_A(T, p)
        B = self.get_B(T, p)
        C = self.get_C(T, p)

        A1 = self.get_pdiff_A_p_T(T, p)
        B1 = self.get_pdiff_B_p_T(T, p)
        C1 = self.get_pdiff_C_p_T(T, p)

        return -(Z**2*A1 + Z*B1 + C1)/(2*A*Z + B + 3*Z**2)

    # Departure function for Cp
    def get_departure_cp(self, T, p, **kwargs):
        """
        Get the departure (difference between real gas and ideal gas) for isobaric specific heat capacity (C_p) [J/mol/K]

        Parameters
        ----------
        T - Temperature [K]
        p - Pressure [Pa]

        Optional parameters
        -------------------
        Z - Compressibility factor [dimensionless]
        This function recalculates compressibility factor if it is not given as an optional parameter

        Returns
        -------
        Departure for isobaric specific heat capacity [J/mol/K]
        """
        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        a = self.get_a(T)

        Z1 = self.get_pdiff_Z_T_p(T, p, Z=Z)
        a1 = self.get_diff_a_T(T)
        a2 = self.get_double_diff_a_T(T)

        D = (self.d**2 - 4*self.e)**0.5
        h = (D + 2*self.R*T*Z/p + self.d)/(-D + 2*self.R*T*Z/p + self.d)
        h1 = -4*D*self.R*p*(T*Z1 + Z)/(2*self.R*T*Z - p*(D - self.d))**2

        return self.R*T*Z1 + self.R*Z - self.R + T*log(h)*a2/D + (T*a1 - a)*h1/D*h

    # Departure function for Cv
    def get_departure_cv(self, T, p, **kwargs):
        """
        Get the departure (difference between real gas and ideal gas) for isochoric specific heat capacity (C_p) [J/mol/K]

        Parameters
        ----------
        T - Temperature [K]
        p - Pressure [Pa]

        Optional parameters
        -------------------
        Z - Compressibility factor [dimensionless]
        This function recalculates compressibility factor if it is not given as an optional parameter

        Returns
        -------
        Departure for isochoric specific heat capacity [J/mol/K]
        """
        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        a2 = self.get_double_diff_a_T(T)

        D = (self.d**2 - 4*self.e)**0.5
        h = (D + 2*self.R*T*Z/p + self.d)/(-D + 2*self.R*T*Z/p + self.d)

        return T/D*a2 * log(h)

    # Isothermal compressibililty
    def get_isothermal_compressibility(self, T, p, **kwargs):
        """
        Get the isothermal compressibility of a real gas

        Parameters
        ----------
        T - Temperature [K]
        p - Pressure [Pa]

        Optional parameters
        -------------------
        Z - Compressibility factor [dimensionless]
        This function recalculates compressibility factor if it is not given as an optional parameter

        Returns
        -------
        Isothermal compressibility [1/Pa]
        """

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        Z1 = self.get_pdiff_Z_p_T(T, p, Z=Z)

        return 1/p - 1/Z * Z1
