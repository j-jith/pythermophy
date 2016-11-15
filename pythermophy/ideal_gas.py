from __future__ import division, print_function

from parent_class import EOS

class IdealGas(EOS):

    def __init__(self, fluid):
        """
        Ideal gas equation of state

        Parameters
        ----------
        fluid - A Fluid class object

        Returns
        -------
        An equation of state object
        """
        super(IdealGas, self).__init__(fluid)

    def get_isothermal_compressibility(self, T, p, **kwargs):
        """
        Get the isothermal compressibility of ideal gas

        Parameters
        ----------
        T - Temperature [K] (only required for consistency with the same function for other classes; not used in computation)
        p - Pressure [Pa]

        Returns
        -------
        Isothermal compressibility [1/Pa]
        """
        return 1/p
