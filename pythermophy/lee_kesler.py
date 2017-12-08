from __future__ import division, print_function
import numpy as np
import scipy.optimize as spo
from math import exp

from .parent_class import EOS

class LeeKesler(EOS):
    """
    Lee-Kesler equation of state.

    :param fluid: a :class:`~pythermophy.fluid.Fluid` instance

    :return: a :class:`~pythrmophy.parent_class.EOS` instance
    """
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
        """
        Lee-Kesler equation of state
        For details see: Chemical Engineering Thermodynamics, Y.V.C. Rao, 1997, University Press (India), pp.79
        https://books.google.co.in/books?id=GjlO9MA9edUC&pg=PA79&dq=lee+kesler+method

        Parameters
        ----------
        fluid - A Fluid class object

        Returns
        -------
        An equation of state object
        """

        super(LeeKesler, self).__init__(fluid)
        self.acentric = fluid.acentric
        self.p_crit = fluid.p_crit # Pa
        self.T_crit = fluid.T_crit # K

    def get_B(self, Tr, fluid):
        """
        Returns :math:`B(T_r) = b_1 - b_2/T_r - b_3/T_r^2 - b_4/T_r^3`.

        :param Tr: reduced temperature (:math:`T_r = T/T_c`)
        :type Tr: float
        :param fluid: specifies if the fluid is *simple* or *reference*
        :type fluid: str

        :return: :math:`B(T_r)`
        :rtype: float
        """
        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        return self.b1[i] - self.b2[i]/Tr - self.b3[i]/Tr**2 - self.b4[i]/Tr**3

    def get_C(self, Tr, fluid):
        """
        Returns :math:`C(T_r) = c_1 - c_2/T_r + c_3/T_r^3`.

        :param Tr: reduced temperature (:math:`T_r = T/T_c`)
        :type Tr: float
        :param fluid: specifies if the fluid is *simple* or *reference*
        :type fluid: str

        :return: :math:`C(T_r)`
        :rtype: float
        """
        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        return self.c1[i] - self.c2[i]/Tr + self.c3[i]/Tr**3

    def get_D(self, Tr, fluid):
        """
        Returns :math:`D(T_r) = d_1 + c_2/T_r`.

        :param Tr: reduced temperature (:math:`T_r = T/T_c`)
        :type Tr: float
        :param fluid: specifies if the fluid is *simple* or *reference*
        :type fluid: str

        :return: :math:`D(T_r)`
        :rtype: float
        """
        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        return self.d1[i] + self.d2[i]/Tr

    def get_reduced_volume(self, Tr, pr, fluid, **kwargs):
        """
        Returns the reduced volume :math:`v_r`.

        This is done by solving the following nonlinear equation:

        .. math::

            \\frac{p_r v_r}{T_r} = 1 + \\frac{B(T_r)}{v_r} + \\frac{C(T_r)}{v_r^2}
            + \\frac{D(T_r)}{v_r^5} + \\frac{c_4}{T_r^3 v_r^2}
            (\\beta + \\frac{\gamma}{v_r^2}) \exp(-\\frac{\gamma}{v_r^2})

        :param Tr: reduced temperature (:math:`T_r = T/T_c`)
        :type Tr: float
        :param pr: reduced pressure (:math:`p_r = p/p_c`)
        :type pr: float
        :param fluid: specifies if the fluid is *simple* or *reference*
        :type fluid: str
        :keyword init_vol: (optional kwarg) initial guess for reduced volume to be used in the
            nonlinear solver
        :type init_vol: float

        :return: reduced volume
        :rtype: float
        """
        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        B = self.get_B(Tr, fluid)
        C = self.get_C(Tr, fluid)
        D = self.get_D(Tr, fluid)

        def objfun(x):
            return -pr*x/Tr + 1 + B/x + C/x**2 + D/x**5 + self.c4[i]/Tr**3/x**2 * (self.beta[i] + self.gamma[i]/x**2) * exp(-self.gamma[i]/x**2)

        # initial gues for nonlinear solver
        if 'init_vol' in kwargs:
            init_vol = kwargs['init_vol']
        else:
            init_vol = Tr/pr
            #if pr<1:
            #    init_vol = 10.
            #else:
            #    init_vol = 0.1
            #if Tr/pr > 1:
            #    init_vol = 10.
            #else:
            #    init_vol = Tr/pr

        #vol, info, ier, mesg = fsolve(objfun, init_vol)
        result = spo.root(objfun, init_vol, method='lm')
        #print(result)

        if len(result.x) == 1:
            return result.x[0]
        else:
            return result.x


    def get_Z(self, T, p, **kwargs):
        """
        Returns the compressibility factor of a real gas.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float
        :keyword init_vol: (optional kwarg) initial guess for
            reduced volumes of ideal and reference fluids to be used in the
            nonlinear equation solver in :func:`get_reduced_volume`. Should be
            a list [v0, vr], where v0 is the inital guess for the simple fluid,
            and vr is the initial guess for the reference fluid.
        :type init_vol: list(float)

        :return: compressibility factor
        :rtype: float
        """

        Tr = T/self.T_crit
        pr = p/self.p_crit

        if 'init_vol' in kwargs:
            init_vol = kwargs['init_vol']
        else:
            init_vol = [Tr/pr, Tr/pr]

        z0 = pr/Tr * self.get_reduced_volume(Tr, pr, 'simple', init_vol=init_vol[0])
        zr = pr/Tr * self.get_reduced_volume(Tr, pr, 'reference', init_vol=init_vol[1])

        # departure term
        z1 = (zr - z0)/self.acentric_r

        z = z0 + self.acentric*z1

        #print('z0 =', z0)
        #print('z1 =', z1)

        return z


    def get_departure_cv_aux(self, Tr, pr, fluid, **kwargs):
        """
        Auxiliary function used by :func:`get_departure_cv`.

        Computes the departure in isochoric specific heat capacity for *simple*
        and *reference* fluids.

        :param Tr: reduced temperature
        :type Tr: float
        :param pr: reduced pressure
        :type pr: float
        :param fluid: specifices if the fluid is *simple* or *reference*
        :type fluid: str
        :keyword Vr: (optional kwarg) reduced volume at (Tr, pr). If not
            provided, :func:`get_reduced_volume` is used to compute it.
        :type Vr: float
        """
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

        E = self.c4[i]/2/Tr**3/self.gamma[i] * ( self.beta[i] + 1 - (self.beta[i] + 1 + self.gamma[i]/vr**2)*exp(-self.gamma[i]/vr**2) )

        d_cv = 2*(self.b3[i] + 3*self.b4[i]/Tr)/Tr**2/vr - 3*self.c3[i]/Tr**3/vr**2 - 6*E

        return d_cv

    def get_departure_cv(self, T, p):
        """
        Returns the departure (difference between real gas and ideal gas) of
        isochoric specific heat capacity :math:`c_v` (J/mol/K)

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: departure of isochoric specific heat capacity (J/mol/K)
        :rtype: float
        """

        Tr = T/self.T_crit
        pr = p/self.p_crit

        cv_0 = self.get_departure_cv_aux(Tr, pr, 'simple')
        cv_r = self.get_departure_cv_aux(Tr, pr, 'reference')

        cv_1 = (cv_r - cv_0)/self.acentric_r

        return self.R * (cv_0 + self.acentric*cv_1)



    def get_pdiff_pr_Tr_Vr(self, Tr, Vr, fluid):
        """
        Returns the partial derivative of :math:`p_r` wrt. :math:`T_r` at
        constant :math:`v_r`.

        :param Tr: reduced temperature
        :type Tr: float
        :param Vr: reduced volume
        :type Vr: float
        :param fluid: specifices if the fluid is *simple* or *reference*
        :type fluid: str

        :rtype: float
        """
        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        # pdiff = 1/Vr * ( 1 + (self.b1[i] + self.b3[i]/Tr**2 + 2*self.b4[i]/Tr**3)/Vr
        #                  + (self.c1[i] - 2*self.c3[i]/Tr**3)/2/Vr**2 + self.d1[i]/5/Vr**5
        #                  - 2*self.c4[i]/Tr**3/Vr**2 * (self.beta[i] + self.gamma[i]/Vr**2) * exp(-self.gamma[i]/Vr**2) )
        pdiff = 1/Vr * ( 1 + (self.b1[i] + self.b3[i]/Tr**2 + 2*self.b4[i]/Tr**3)/Vr
                         + (self.c1[i] - 2*self.c3[i]/Tr**3)/Vr**2 + self.d1[i]/Vr**5
                         - 2*self.c4[i]/Tr**3/Vr**2 * (self.beta[i] + self.gamma[i]/Vr**2) * exp(-self.gamma[i]/Vr**2) )

        return pdiff

    def get_pdiff_pr_Vr_Tr(self, Tr, Vr, fluid):
        """
        Returns the partial derivative of :math:`p_r` wrt. :math:`V_r` at
        constant :math:`T_r`.

        :param Tr: reduced temperature
        :type Tr: float
        :param Vr: reduced volume
        :type Vr: float
        :param fluid: specifices if the fluid is *simple* or *reference*
        :type fluid: str

        :rtype: float
        """
        if fluid=='reference':
            i = 1
        elif fluid=='simple':
            i = 0
        else:
            return None

        B = self.get_B(Tr, fluid)
        C = self.get_C(Tr, fluid)
        D = self.get_D(Tr, fluid)

        pdiff = -Tr/Vr**2 * ( 1 + 2*B/Vr + 3*C/Vr**2 + 6*D/Vr**5
                              + self.c4[i]/Tr**3/Vr**2
                              * (3*self.beta[i] + (5-2*self.beta[i]-2*self.gamma[i]/Vr**2)*self.gamma[i]/Vr**2)
                              * exp(-self.gamma[i]/Vr**2) )

        return pdiff

    # Isothermal compressibililty
    def get_isothermal_compressibility(self, T, p):
        """
        Returns the isothermal compressibility of a real gas.

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: isothermal compressibility (1/Pa)
        :rtype: float
        """

        Tr = T/self.T_crit
        pr = p/self.p_crit

        fluids = ['simple', 'reference']

        beta_T = []


        for i in [0,1]:
            vr = self.get_reduced_volume(Tr, pr, fluids[i])
            dv_dp = 1./self.get_pdiff_pr_Vr_Tr(Tr, vr, fluids[i])
            beta_T.append(-dv_dp/vr / self.p_crit)

        beta1 = (beta_T[1] - beta_T[0])/self.acentric_r
        return beta_T[0] + self.acentric*beta1


    def get_departure_cp_aux(self, Tr, pr, fluid):
        """
        Auxiliary function used by :func:`get_departure_cp`.

        Computes the departure in isobaric specific heat capacity for *simple*
        and *reference* fluids.

        :param Tr: reduced temperature
        :type Tr: float
        :param pr: reduced pressure
        :type pr: float
        :param fluid: specifices if the fluid is *simple* or *reference*
        :type fluid: str
        """

        vr = self.get_reduced_volume(Tr, pr, fluid)

        cv = self.get_departure_cv_aux(Tr, pr, fluid, Vr=vr)

        dp_dT = self.get_pdiff_pr_Tr_Vr(Tr, vr, fluid)
        dp_dv = self.get_pdiff_pr_Vr_Tr(Tr, vr, fluid)

        #print(dp_dv)
        #print(dp_dT)

        cp = cv - 1 - Tr * dp_dT**2 / dp_dv

        return cp

    def get_departure_cp(self, T, p):
        """
        Returns the departure (difference between real gas and ideal gas) of
        isobaric specific heat capacity :math:`c_p` (J/mol/K)

        :param T: temperature (K)
        :type T: float
        :param p: pressure (Pa)
        :type p: float

        :return: departure of isobaric specific heat capacity (J/mol/K)
        :rtype: float
        """

        Tr = T/self.T_crit
        pr = p/self.p_crit

        cp_0 = self.get_departure_cp_aux(Tr, pr, 'simple')
        cp_r = self.get_departure_cp_aux(Tr, pr, 'reference')

        cp_1 = (cp_r - cp_0)/self.acentric_r

        #print('cp_0 - cp* =', cp_0)
        #print('cp_1 - cp* =', cp_1)

        return self.R * (cp_0 + self.acentric*cp_1)


        # vr = self.get_reduced_volume(Tr, pr, 'simple')
        # dp_dT_0 = self.get_pdiff_pr_Tr_Vr(Tr, vr, 'simple')
        # dp_dv_0 = self.get_pdiff_pr_Vr_Tr(Tr, vr, 'simple')
        # cv_0 = self.get_departure_cv_aux(Tr, pr, 'simple', Vr=vr)

        # vr = self.get_reduced_volume(Tr, pr, 'reference')
        # dp_dT_r = self.get_pdiff_pr_Tr_Vr(Tr, vr, 'reference')
        # dp_dv_r = self.get_pdiff_pr_Vr_Tr(Tr, vr, 'reference')
        # cv_r = self.get_departure_cv_aux(Tr, pr, 'reference', Vr=vr)

        # #dp_dT = dp_dT_0 + self.acentric/self.acentric_r*(dp_dT_r - dp_dT_0)
        # #dp_dv = dp_dv_0 + self.acentric/self.acentric_r*(dp_dv_r - dp_dv_0)
        # #cv = cv_0 + self.acentric/self.acentric_r*(cv_r - cv_0)

        # cp_0 = cv_0 - 1 - Tr*dp_dT_0**2/dp_dv_0
        # cp_r = cv_r - 1 - Tr*dp_dT_r**2/dp_dv_r

        # cp_1 = (cp_r - cp_0)/self.acentric_r

        # print('cp_0 - cp* =', cp_0)
        # print('cp_1 - cp* =', cp_1)

        # cp = self.R*(cp_0 + self.acentric*cp_1)

        # return cp
