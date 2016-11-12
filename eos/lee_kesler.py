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

    def __init__(self, fluid):

        super(LeeKesler, self).__init__(fluid)
        self.acentric = fluid.acentric
        self.p_crit = fluid.p_crit # Pa
        self.T_crit = fluid.T_crit # K

    def get_reduced_volume(self, Tr, pr, fluid='simple'):

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

        if len(result.x) == 1:
            return result.x[0]
        else:
            return result.x

    def get_Z(self, T, p):

        Tr = T/self.T_crit
        pr = p/self.p_crit

        z0 = pr/Tr * self.get_reduced_volume(Tr, pr, 'simple')
        zr = pr/Tr * self.get_reduced_volume(Tr, pr, 'reference')

        # departure term
        z1 = (zr - z0)/self.acentric_r

        z = z0 + self.acentric*z1

        #print('z0 = ', z0)
        #print('zr = ', zr)
        #print('z1 = ', z1)
        #print('z = ', z)

        return z

    def get_reduced_departure_energy(self, Tr, pr, fluid='simple', **kwargs):
        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        if 'Vr' in kwargs:
            vr = kwargs['Vr']
        else:
            vr = self.get_reduced_volume(Tr, pr, fluid)

        du = Tr*(-1/Tr/vr * (self.b2[i] + 2*self.b3[i]/Tr + self.b4[i]/Tr**2) \
                -1/2/Tr/vr**2 * (self.c2[i] - 3*self.c3[i]/Tr**2) \
                + self.d2[i]/5/Tr/vr**5 \
                + 3*self.c4[i]/2/Tr**3/self.gamma[i] * \
                (self.beta[i] + 1 - (self.beta[i] + 1 + self.gamma[i]/vr**2)*exp(-self.gamma[i]/vr**2)))

        return du


    def get_departure_energy(self, T, p, **kwargs):

        Tr = T/self.T_crit
        pr = p/self.p_crit

        if 'Vr' in kwargs:
            vr = kwargs['Vr']
            du0 = self.R*self.T_crit * self.get_reduced_departure_energy(Tr, pr, 'simple', Vr=vr[0])
            dur = self.R*self.T_crit * self.get_reduced_departure_energy(Tr, pr, 'reference', Vr=vr[1])
        else:
            du0 = self.R*self.T_crit * self.get_reduced_departure_energy(Tr, pr, 'simple')
            dur = self.R*self.T_crit * self.get_reduced_departure_energy(Tr, pr, 'reference')

        # departure term
        du1 = (dur - du0)/self.acentric_r

        return du0 + self.acentric*du1


    def get_departure_enthalpy(self, T, p, **kwargs):

        u = self.get_departure_energy(T, p)

        if 'Z' in kwargs:
            z = kwargs['Z']
        else:
            z = self.get_Z(T, p)

        return self.R*T*(z-1) + u

    def get_departure_cp(self, T, p, step=1e-3, **kwargs):

        h2 = self.get_departure_enthalpy(T+step, p, **kwargs)
        h1 = self.get_departure_enthalpy(T-step, p, **kwargs)

        return (h2 - h1)/2/step

    def get_departure_cv(self, T, p, step=1e-3, **kwargs):
        Tr = T/self.T_crit
        pr = p/self.p_crit

        vr0 = self.get_reduced_volume(Tr, pr, 'simple')
        vrr = self.get_reduced_volume(Tr, pr, 'reference')

        u2 = self.get_departure_energy(T+step, p, Vr=[vr0, vrr])
        u1 = self.get_departure_energy(T-step, p, Vr=[vr0, vrr])

        return (u2 - u1)/2/step

    def get_pdiff_Z_p_T(self, T, p, step=1e-3):
        z2 = self.get_Z(T, p+step)
        z1 = self.get_Z(T, p-step)
        return (z2 - z1)/2/step

    # Isothermal compressibililty
    def get_isothermal_compressibility(self, T, p, **kwargs):

        if 'Z' in kwargs:
            Z = kwargs['Z']
        else:
            Z = self.get_Z(T, p)

        Z1 = self.get_pdiff_Z_p_T(T, p)

        return 1/p - 1/Z * Z1
