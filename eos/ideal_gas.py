from __future__ import division, print_function

from parent_class import EOS

class IdealGas(EOS):

    def __init__(self, M, fluid):
        super(IdealGas, self).__init__(M, fluid)
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

    def get_isothermal_compressibility(self, T, p, **kwargs):
        return 1/p
