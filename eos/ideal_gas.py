from __future__ import division, print_function

from parent_class import EOS

class IdealGas(EOS):

    def __init__(self, fluid):
        super(IdealGas, self).__init__(fluid)

    def get_isothermal_compressibility(self, T, p, **kwargs):
        return 1/p
