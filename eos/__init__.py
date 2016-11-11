#from parent_class import EOS
from ideal_gas import IdealGas

from cubic_eos import RedlichKwong
from cubic_eos import SoaveRedlichKwong
from cubic_eos import PengRobinson

from lee_kesler import LeeKesler

# Shorter aliases for class names
class RK(RedlichKwong):
    pass

class SRK(SoaveRedlichKwong):
    pass

class PR(PengRobinson):
    pass

class LK(LeeKesler):
    pass

