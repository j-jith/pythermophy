from parent_class import EOS
from ideal_gas import IdealGas
from redlich_kwong import RedlichKwong
from soave_redlich_kwong import SoaveRedlichKwong
from peng_robinson import PengRobinson
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

