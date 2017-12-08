from __future__ import division, print_function

class EOS(object):
    """
    The parent equation of state (EOS) class. All the other EOS classes are
    derived from this.

    :param fluid: a :class:`~pythermophy.fluid.Fluid` instance

    :return: a :class:`EOS` instance
    """

    # Universal gas constant
    R = 8.3144598 # J/K/mol

    def __init__(self, fluid):
        self.fluid = fluid
        self.molar_mass = fluid.molar_mass
        self.R_sp = self.R/self.molar_mass
        self.ideal_cp_coeffs = fluid.ideal_cp_coeffs

    def get_Z(self, T, p):
        """
        Returns the compressibility factor of the fluid.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: compressibility factor
        :rtype: float
        """
        return 1.

    def get_rho(self, T, p, **kwargs):
        """
        Returns the density of the fluid.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float
        """
        if 'Z' in kwargs:
            z = kwargs['Z']
        else:
            z = self.get_Z(T, p)
        #print(z)
        return p/z/self.R_sp/T

    def get_ideal_cp(self, T):
        """
        Returns the ideal gas isobaric specific heat capacity (:math:`c_p^0`) in J/mol/K.
        See http://www.wiley.com/college/moran/CL_0471465704_S/user/tables/TABLE3S/table3sframe.html

        :param T: temperature (K)
        :type T: float

        :return: ideal gas isobaric specific heat capacity (J/mol/K)
        :rtype: float
        """

        cp = 0.
        for i, c_i in enumerate(self.ideal_cp_coeffs):
            cp = cp + c_i*T**i

        # try:
        #     if self.fluid=='CO2':

        #         c = [22.26, 5.981e-2, -3.501e-5, 7.469e-9]
        #         for i, c_i in enumerate(c):
        #             cp = cp + c_i*T**i

        #     else:
        #         raise ValueError('Ideal gas Cp for fluid "{}" not implemented'.format(self.fluid))
        # except ValueError as err:
        #     print(err)

        return cp

    def get_ideal_cv(self, T):
        """
        Returns the ideal gas isochoric specific heat capacity (:math:`c_v^0`)
        in J/mol/K.

        :param T: temperature (K)
        :type T: float

        :return: ideal gas isochoric specific heat capacity (J/mol/K)
        :rtype: float
        """
        return self.get_ideal_cp(T) - self.R

    def get_departure_cp(self, T, p):
        """
        Returns the departure (difference between real gas and ideal gas) of
        isobaric specific heat capacity (:math:`c_p`) in J/mol/K.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: departure for isobaric specific heat capacity (J/mol/K)
        :rtype: float
        """
        return 0.

    def get_departure_cv(self, T, p):
        """
        Returns the departure (difference between real gas and ideal gas) of
        isochoric specific heat capacity (:math:`c_v`) in J/mol/K.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: departure for isochoric specific heat capacity (J/mol/K)
        :rtype: float
        """
        return 0.

    def get_cp(self, T, p, **kwargs):
        """
        Returns the real gas isobaric specific heat capacity (:math:`c_p`) in
        J/mol/K.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: isobaric specific heat capacity (J/mol/K)
        :rtype: float
        """
        return self.get_ideal_cp(T) + self.get_departure_cp(T, p, **kwargs)

    def get_cv(self, T, p, **kwargs):
        """
        Returns the real gas isochoric specific heat capacity (:math:`c_v`) in
        J/mol/K.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: isochoric specific heat capacity (J/mol/K)
        :rtype: float
        """
        return self.get_ideal_cv(T) + self.get_departure_cv(T, p, **kwargs)

    def get_adiabatic_index(self, T, p, **kwargs):
        """
        Returns the adiabatic index (ratio between isobaric and isochoric heat
        capacities) for a real gas.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: adiabatic index (dimensionless)
        :rtype: float
        """
        return self.get_cp(T, p, **kwargs)/self.get_cv(T, p, **kwargs)

    def get_speed_of_sound(self, T, p, **kwargs):
        """
        Returns the speed of sound in a real gas.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: speed of sound (m/s)
        :rtype: float
        """
        if 'Z' in kwargs:
            z = kwargs['Z']
        else:
            z = self.get_Z(T, p)

        gamma = self.get_adiabatic_index(T, p, **kwargs)
        beta = self.get_isothermal_compressibility(T, p, **kwargs)
        rho = self.get_rho(T, p, **kwargs)

        c = (gamma/beta/rho)**0.5

        return c
