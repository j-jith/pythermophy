from __future__ import print_function, division
import numpy as np
from math import log

from .parent_class import EOS

class CubicEOS(EOS):
    """
    Parent class from which other cubic equations of state are derived.

    The equation of state is expressed as:

    .. math::

        p = \\frac{RT}{v - b + c} - \\frac{a(T)}{v^2 + dv + e}

    where :math:`(p, v, T)` are the pressure, specific volume, and temperature
    respectively. :math:`R` is the specific gas constant, and :math:`(a(T), b,
    c, d, e)` are the coefficients of the cubic equation of state.

    :param b: coefficient :math:`b`
    :type b: float
    :param c: coefficient :math:`c`
    :type c: float
    :param d: coefficient :math:`d`
    :type d: float
    :param e: coefficient :math:`e`
    :type e: float
    :param fluid: a :class:`~pythermophy.fluid.Fluid` instance

    :return: a :class:`~pythermophy.parent_class.EOS` instance
    """

    def __init__(self, b, c, d, e, fluid):

        self.b = b; self.c = c
        self.d = d; self.e = e

        super(CubicEOS, self).__init__(fluid)

    def get_a(self, T):
        """
        Returns the temperature dependent coefficient :math:`a(T)`.

        Function should be redefined in the child class.
        """
        pass

    def get_diff_a_T(self, T):
        """
        Returns the derivative of coefficient :math:`a(T)` wrt. temperature :math:`T`.

        Function should be redefined in the child class.
        """
        pass

    def get_double_diff_a_T(self, T):
        """
        Returns the second derivative of coefficient :math:`a(T)` wrt. temperature :math:`T`.

        Function should be redefined in the child class.
        """
        pass

    # Coefficients of the cubic polynomial for Z: A, B, C
    # Z^3 + A*Z^2 + B*Z + C = 0
    def get_A(self, T, p):
        """
        Returns the coefficient :math:`A` of the cubic polynomial of
        compressiblity factor given by :math:`Z^3 + A Z^2 + B Z + C = 0`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: :math:`A`
        :rtype: float
        """
        return (-self.R*T - self.b*p + self.c*p + self.d*p)/(self.R*T)

    def get_B(self, T, p):
        """
        Returns the coefficient :math:`B` of the cubic polynomial of
        compressiblity factor given by :math:`Z^3 + A Z^2 + B Z + C = 0`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: :math:`B`
        :rtype: float
        """
        a = self.get_a(T)
        return (-self.R*T*self.d*p + a*p - self.b*self.d*p**2 + self.c*self.d*p**2 + self.e*p**2)/(self.R**2*T**2)

    def get_C(self, T, p):
        """
        Returns the coefficient :math:`C` of the cubic polynomial of
        compressiblity factor given by :math:`Z^3 + A Z^2 + B Z + C = 0`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: :math:`C`
        :rtype: float
        """
        a = self.get_a(T)
        return (-self.R*T*self.e*p**2 - a*self.b*p**2 + a*self.c*p**2 - self.b*self.e*p**3 + self.c*self.e*p**3)/(self.R**3*T**3)

    # Partial derivatives of A, B and C wrt. T at constant p
    def get_pdiff_A_T_p(self, T, p):
        """
        Returns the partial derivative of coefficient :math:`A` wrt. :math:`T`
        at constant :math:`p`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :rtype: float
        """
        return p*(self.b - self.c - self.d)/(self.R*T**2)

    def get_pdiff_B_T_p(self, T, p):
        """
        Returns the partial derivative of coefficient :math:`B` wrt. :math:`T`
        at constant :math:`p`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :rtype: float
        """
        a = self.get_a(T)
        dadT = self.get_diff_a_T(T)
        return p*(self.R*T*self.d + T*dadT + 2*self.b*self.d*p - 2*self.c*self.d*p - 2*self.e*p - 2*a)/(self.R**2*T**3)

    def get_pdiff_C_T_p(self, T, p):
        """
        Returns the partial derivative of coefficient :math:`C` wrt. :math:`T`
        at constant :math:`p`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :rtype: float
        """
        a = self.get_a(T)
        dadT = self.get_diff_a_T(T)
        return p**2*(2*self.R*T*self.e + T*(-self.b + self.c)*dadT + 3*self.b*self.e*p + 3*self.b*a - 3*self.c*self.e*p - 3*self.c*a)/(self.R**3*T**4)

    # Partial derivatives of A, B and C wrt. p at constant T
    def get_pdiff_A_p_T(self, T, p):
        """
        Returns the partial derivative of coefficient :math:`A` wrt. :math:`p`
        at constant :math:`T`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :rtype: float
        """
        return (-self.b + self.c + self.d)/(self.R*T)

    def get_pdiff_B_p_T(self, T, p):
        """
        Returns the partial derivative of coefficient :math:`B` wrt. :math:`p`
        at constant :math:`T`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :rtype: float
        """
        a = self.get_a(T)
        return (-self.R*T*self.d - 2*self.b*self.d*p + 2*self.c*self.d*p + 2*self.e*p + a)/(self.R**2*T**2)

    def get_pdiff_C_p_T(self, T, p):
        """
        Returns the partial derivative of coefficient :math:`C` wrt. :math:`p`
        at constant :math:`T`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :rtype: float
        """
        a = self.get_a(T)
        return p*(-2*self.R*T*self.e - 3*self.b*self.e*p - 2*self.b*a + 3*self.c*self.e*p + 2*self.c*a)/(self.R**3*T**3)


    # Solve the cubic equation for compressiblilty factor Z
    def get_Z(self, T, p):
        """
        Returns the compressibility factor of a real gas.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: compressibility factor
        :rtype: float
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
        """
        Returns the partial derivative of :math:`Z` wrt. :math:`T` at constant
        :math:`p`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :keyword Z: (optional kwarg) compressibility factor at (T, p). Computed
            by :func:`get_Z` if not provided.
        :type Z: float

        :rtype: float
        """
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
        """
        Returns the partial derivative of :math:`Z` wrt. :math:`p` at constant
        :math:`T`.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :keyword Z: (optional kwarg) compressibility factor at (T, p). Computed
            by :func:`get_Z` if not provided.
        :type Z: float

        :rtype: float
        """
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
        Returns the departure (difference between real gas and ideal gas) of
        isobaric specific heat capacity :math:`c_p` (J/mol/K)

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float
        :keyword Z: (optional kwarg) compressibility factor at (T, p). Computed
            by :func:`get_Z` if not provided.
        :type Z: float

        :return: departure of isobaric specific heat capacity (J/mol/K)
        :rtype: float
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
        Returns the departure (difference between real gas and ideal gas) for
        isochoric specific heat capacity (:math:`c_v`) in J/mol/K.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float
        :keyword Z: (optional kwarg) compressibility factor at (T, p). Computed
            by :func:`get_Z` if not provided.
        :type Z: float

        :return: departure of isochoric specific heat capacity (J/mol/K)
        :rtype: float
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
        Returns the isothermal compressibility of a real gas.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: isothermal compressibility (1/Pa)
        :rtype: float
        """

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        Z1 = self.get_pdiff_Z_p_T(T, p, Z=Z)

        return 1/p - 1/Z * Z1
