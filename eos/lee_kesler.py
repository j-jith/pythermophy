from __future__ import division, print_function
import numpy as np
import scipy.optimize as spo
from math import exp, log

from parent_class import EOS

class LeeKesler(EOS):
    '''https://books.google.co.in/books?id=GjlO9MA9edUC&pg=PA79&dq=lee+kesler+method
       Chemical Engineering Thermodynamics, Y.V.C. Rao, 1997, University Press (India), pp.79
    '''

    b1 = [0.1181193, 0.2026579]
    b2 = [0.265728, 0.331511]
    b3 = [0.154790, 0.027655]
    b4 = [0.030323, 0.203488]

    c1 = [0.0236744, 0.0313385]
    c2 = [0.0186984, 0.0503618]
    c3 = [0.0, 0.016901]
    c4 = [0.042724, 0.041577]

    d1 = [0.155488e-4, 0.48736e-4]
    d2 = [0.623689e-4, 0.0740336e-4]

    beta = [0.65392, 1.226]
    gamma = [0.060167, 0.03754]

    acentric_r = 0.3978

    def __init__(self, Tc, pc, M, omega, fluid):

        super(LeeKesler, self).__init__(M, fluid)
        self.acentric = omega
        self.p_crit = pc # Pa
        self.T_crit = Tc # K
        self.molar_mass = M # kg/mol
        self.R_sp = self.R/M

    def get_reduced_vol(self, Tr, pr, fluid='simple'):

        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        B = self.b1[i] - self.b2[i]/Tr - self.b3[i]/Tr**2 - self.b4[i]/Tr**3
        C = self.c1[i] - self.c2[i]/Tr + self.c3[i]/Tr**3
        D = self.d1[i] + self.d2[i]/Tr

        def objfun(x):
            return -pr*x/Tr + 1 + B/x + C/x**2 + D/x**5 + self.c4[i]/Tr**3/x**2 * (self.beta[i] + self.gamma[i]/x**2) * exp(-self.gamma[i]/x**2)

        # initial gues for nonlinear solver
        if pr<1:
            init_vol = 10.
        else:
            init_vol = 0.1

        #vol, info, ier, mesg = fsolve(objfun, init_vol)
        result = spo.root(objfun, init_vol, method='lm')
        #print(result)

        return result.x

    def get_Z(self, T, p):

        Tr = T/self.T_crit
        pr = p/self.p_crit

        z0 = pr/Tr * self.get_reduced_vol(Tr, pr, 'simple')
        zr = pr/Tr * self.get_reduced_vol(Tr, pr, 'reference')

        # departure term
        z1 = (zr - z0)/self.acentric_r

        z = z0 + self.acentric*z1

        #print('z0 = ', z0)
        #print('zr = ', zr)
        #print('z1 = ', z1)
        #print('z = ', z)

        return z

