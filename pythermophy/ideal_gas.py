from __future__ import division, print_function

from .parent_class import EOS

class IdealGas(EOS):
    """
    Ideal gas equation of state.

    :param fluid: a :class:`~pythermophy.fluid.Fluid` instance

    :return: an equation of state (:class:`~pythermophy.parent_class.EOS`) instance
    """

    def __init__(self, fluid):
        super(IdealGas, self).__init__(fluid)

    def get_isothermal_compressibility(self, T, p, **kwargs):
        """
        Returns the isothermal compressibility of an ideal gas

        :param T: temperature (K) [only required for consistency with the
            same function for other classes; not used in computation]
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: isothermal compressibility (1/Pa)
        :rtype: float
        """
        return 1/p
