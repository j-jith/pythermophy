Theory
######

.. default-role:: math

Equation of state (EOS) of a real fluid can be expressed in a general form as

.. math::

   pv = ZRT,

where `p` is the pressure, `v` is the specific volume, `T` is the temperature,
`Z` is the compressibility factor, and `R = R_0 / M` is the specific gas
constant (`R_0` is the universal gas constant and `M` is the molar mass of the
fluid). Various EOS differ based on how they determine `Z`.

The density of the fluid can be easily obtained from the specific volume as
`\rho=M/v`. The speed of sound, on the other hand requires knowledge of the
heat capacities, and isothermal compressibility of the fluid. Tee expression
for speed of sound of a real fluid is

.. math::

    c^2 = \frac{1}{\rho} \frac{c_p}{c_v} \frac{1}{\chi_T},

where `c_p` and `c_v` are the isobaric and isochoric specific heat capacities,
and `chi_t` is the isothermal compressibility.

Ideal Gas EOS
=============

The ideal gas (IG) EOS is obtained by assuming `Z=1`. This EOS is quite
satisfactory at lower pressures. However, at higher pressures it is unreliable
as it does not take into account the finite volume of the fluid molecules or
the interactions between them.

For an ideal gas, the heat capacities can be expressed as polynomials of
temperature. In the present work, the isobaric specific heat capacity is
assumed to be modelled by the following polynomial.

.. math::

    c_p = a + bT + cT^2 + dT^3.

The coefficients `(a, b, c, d)` for various fluids can be found in the tables
provided by ThermoNet_ (reproduced here in :doc:`table_cp`). The ideal gas
isochoric specific heat capacity can then be obtained by the relation `c_p -c_v
= R`.

The isothermal compressibility is of an ideal gas is given by

.. math::

    \chi_T = -\frac{1}{v} \left. \frac{\partial v}{\partial p} \right|_T = \frac{1}{p}

Cubic EOS
=========

The cubic EOS are an improvement on the ideal gas model in the sense that they
try to model the finite volume of the fluid molecules as well as the
interactions between them. The present work considers three different cubic EOS
-- Redlich-Kwong (RK) [redlich1949]_, Soave-Redlich-Kwong (SRK) [soave1972]_,
and Peng-Robinson(PR) [peng1976]_. A general expression for cubic EOS is given
by

.. math::

    p = \frac{RT}{v-b+c} - \frac{a(T)}{v^2+dv+e},

where `(a(T), b, c, d, e)` are the coefficients of the EOS. Substituting the
general EOS `pv = ZRT` into the above equation results in a cubic equation in
terms of `Z` given by

.. math::

    Z^3 + A Z^2 + B Z + C = 0.

The coefficients `A, B, C` are expressed as

.. math::
        A(T, p) = \frac{-RT - bp + cp + dp}{RT} \\
        B(T, p) = \frac{-RTdp + a(T)p - bdp^2 + cdp^2 + ep^2}{R^2 T^2} \\
        C(T, p) = \frac{-RTep^2 - a(T)bp^2 + a(T)cp^2 - bep^3 + cep^3}{R^3
        T^3}.

For the cubic EOS considered here, these coefficients are as follows.

- **Redlich-Kwong (RK)**:

  .. math::

      a(T) = \frac{0.4278 R^2 T_c^2}{p_c} \sqrt{\frac{T_c}{T}} \\
      b = d = \frac{0.08664 R T_c}{p_c} \\
      c = e = 0.

- **Soave-Redlich-Kwong (SRK)**

  .. math::

      a(T) = \frac{0.4278 R^2 T_c^2}{p_c} \left[1 + \kappa \left(1 - \sqrt{\frac{T}{T_c}} \right) \right]^2 \\
      \kappa = 0.48508 + 1.55171 \omega - 0.15613 \omega^2 \\
      b = d = \frac{0.08664 R T_c}{p_c} \\
      c = e = 0.

- **Peng-Robinson (PR)**

  .. math::

      a(T) = \frac{0.45724 R^2 T_c^2}{p_c} \left[1 + \kappa \left(1 - \sqrt{\frac{T}{T_c}} \right) \right]^2 \\
      \kappa = 0.37464 + 1.54226 \omega - 0.26992 \omega^2 \\
      b = \frac{0.07780 R T_c}{p_c} \\
      c = 0 \qquad d = 2b \qquad e = -b^2.

Here, `(T_c, p_c)` are the critical temperature and pressure of the fluid, and
`\omega` is its acentric factor.

The heat capacities can be obtained from the enthalpy and internal energy of
the fluid. The isobaric specific heat capacity `c_p` is given by `c_p =
\left. \frac{\partial H}{\partial T} \right|_p`, where `H` is the specific
enthalpy. Similarly, the isochoric specific heat capacity `c_v` is given by
`c_v = \left. \frac{\partial U}{\partial T} \right|_v`, where `U` is the
specific internal energy. The expressions for `H` and `U` can be found in
[assael1996]_, and are not reproduced here for the sake of brevity. After
incorporating these expressions, the heat capacities can be written as

.. math::

    c_p = c_p^0 %
    + R T Z_p' %
    + R (Z-1) %
    + T \frac{\log(h)}{\Delta} \frac{\mathrm d^2 a}{\mathrm d T^2} %
    + \frac{h g}{\Delta} (T \frac{\mathrm d a}{\mathrm d T} - a) %
    \\
    c_v = c_v^0 %
    + \frac{T}{\Delta} \frac{\mathrm d^2 a}{\mathrm d T^2} \log(h),

where `Z=pv/RT`, and `c_p^0` and `c_v^0` are the ideal gas specific heat
capacities; and

.. math::

    \Delta = \sqrt{d^2 - 4 e} \qquad
    Z_p' = \left. \frac{\partial Z}{\partial T} \right|_p \\
    h = \frac{\Delta + 2 R T Z/p + d}{-\Delta + 2 R T Z/p + d} \qquad
    g = \frac{-4 \Delta R p (T Z_p' + Z)}{(2 R T Z - p(\Delta - d))^2}.

Lastly, the isothermal compressibility of the fluid is obtained using the
following relation.

.. math::

    \chi_T = -\frac{1}{v} \left. \frac{\partial v}{\partial p} \right|_T
        = \frac{1}{p} - \frac{1}{Z} \left. \frac{\partial Z}{\partial p} \right|_T.

Lee-Kesler EOS
==============

The Lee-Kesler (LK) EOS [lee1975]_ is based on the correlation developed by
Pitzer and co. according to which `Z` of a fluid can be written as

.. math::

    Z(T_r, p_r) = Z_0(T_r, p_r) + \omega Z_1 (T_r, p_r),

where `Z_0` is the compressibility factor of a simple fluid whose molecules are
spherical, and `Z_1` is the deviation of compressibility factor. `\omega` is
the acentric factor of the fluid -- a measure of the non-spherical nature of
the molecules. Also, `T_r=T/T_c` and `p_r=p/p_c` are the reduced temperature
and pressure. The deviation `Z_1` is expressed as a linear function of `Z_0`
and the compressibility factor `Z_2` of a heavy non-spherical reference fluid
(with acentric factor `\omega_2`) in the following manner.

.. math::

    Z_1(T_r, p_r) = \frac{Z_2(T_r, p_r) - Z_0(T_r, p_r)}{\omega_2}.

`Z_0` and `Z_2` are obtained by solving the following nonlinear equation.

.. math::

    Z_{0,2} = \frac{p_r v_r}{T_r} = 1 + \frac{B(T_r)}{v_r} + \frac{C(T_r)}{v_r^2}
    + \frac{D(T_r)}{v_r^5}
    + \frac{c_4}{T_r^3 v_r^2} \left( \beta +
    \frac{\gamma}{v_r^2} \right) \exp{\left(\frac{-\gamma}{v_r^2}\right)},

where `v_r` is called the reduced volume, and

.. math::

    B(T_r) = b_1 - \frac{b_2}{T_r} - \frac{b_3}{T_r^2} - \frac{b_4}{T_r^3} \\
    C(T_r) = c_1 - \frac{c_2}{T_r} + {c_3}{T_r^3} \\
    D(T_r) = d_1 + \frac{d_2}{T_r}.

The constant `b_i, c_i, d_i, \beta` and `\gamma` are different for the simple
and reference fluids, and are listed in the table below.

+-------------------+--------------+-----------------+-------------------+--------------+-----------------+
| Constant          | Simple fluid | Reference fluid | Constant          | Simple fluid | Reference fluid |
+===================+==============+=================+===================+==============+=================+
| `b_1`             | 0.1181193    | 0.2026579       | `b_2`             | 0.265728     | 0.331511        |
+-------------------+--------------+-----------------+-------------------+--------------+-----------------+
| `b_3`             | 0.154790     | 0.027655        | `b_4`             | 0.030323     | 0.203488        |
+-------------------+--------------+-----------------+-------------------+--------------+-----------------+
| `c_1`             | 0.0236744    | 0.0313385       | `c_2`             | 0.0186984    | 0.0503618       |
+-------------------+--------------+-----------------+-------------------+--------------+-----------------+
| `c_3`             | 0.0          | 0.016901        | `c_4`             | 0.042724     | 0.041577        |
+-------------------+--------------+-----------------+-------------------+--------------+-----------------+
| `d_1 \times 10^4` | 0.155488$    | 0.48736         | `d_2 \times 10^4` | 0.623689$    | 0.0740336       |
+-------------------+--------------+-----------------+-------------------+--------------+-----------------+
| `\beta`           | 0.65392      | 1.226           | `\gamma`          | 0.060167     | 0.03754         |
+-------------------+--------------+-----------------+-------------------+--------------+-----------------+

The departure (difference between real and ideal fluid) in `c_p` and `c_v` for
the simple and reference fluids is given by [assael1996]_

.. math::

    \frac{\left. c_v^r \right|_{0, 2}}{R} = 2 \frac{(b_3 + 3 b_4/T_r)}{T_r v_r} - 3
        \frac{c_3}{T_r^3 v_r^2} - 6 E \\
        \frac{\left. c_p^r \right|_{0,2}}{R} = \frac{c_v^r|_i}{R} - 1 - T_r
            \left. \left( \left. \frac{\partial p_r}{\partial T_r} \right|_{v_r} \right)^2 \right/ 
                    \left( \left. \frac{\partial p_r}{\partial v_r} \right|_{T_r} \right),

where

.. math::

    E = \frac{c_4}{2 T_r^3 \gamma}  \left( \beta + 1 - \left(\beta + 1 +
    \frac{\gamma}{v_r^2}\right) \exp\left(\frac{-\gamma}{v_r^2}\right) \right).

Then, the real gas heat capacities are given by

.. math::

    c_p = c_p^0 + c_p^r|_0 + \omega \frac{(c_p^r|_2 - c_p^r|_0)}{\omega_2} \\
    c_v = c_v^0 + c_v^r|_0 + \omega \frac{(c_v^r|_2 - c_v^r|_0)}{\omega_2},

where `c_p^0` and `c_v^0` are the ideal gas heat capacities.

The isothermal compressibility of the simple and reference fluids is given by

.. math::

    \chi_T|_{0,2} = -\frac{1}{v} \left. \frac{\partial v}{\partial p} \right|_T
        = \frac{1}{p} - \frac{1}{Z_{0,2}} \left. \frac{\partial Z_{0,2}}{\partial p} \right|_T.

The isothermal compressibility of the real fluid can then be obtained by

.. math::

    \chi_T = \chi_T|_0 + \omega \frac{(\chi_T|_2 - \chi_T|_0)}{\omega_2}.


References
==========

.. [redlich1949] Redlich O, Kwong JNS. On the Thermodynamics of Solutions. V.
   An Equation of State. Fugacities of Gaseous Solutions. Chemical Reviews
   1949;44:233–44. doi:10.1021/cr60137a013.

.. [soave1972] Soave G. Equilibrium constants from a modified Redlich-Kwong
   equation of state. Chemical Engineering Science 1972;27:1197–203.
   doi:10.1016/0009-2509(72)80096-4.

.. [peng1976] Peng D-Y, Robinson DB. A New Two-Constant Equation of State.
   Industrial & Engineering Chemistry Fundamentals 1976;15:59–64.
   doi:10.1021/i160057a011.

.. [assael1996] Assael MJ, Trusler JPM, Tsolakis TF. Thermophysical Properties
   of Fluids: An Introduction to Their Prediction. World Scientific; 1996.

.. [lee1975] Lee BI, Kesler MG. A generalized thermodynamic correlation based
   on three-parameter corresponding states. AIChE J 1975;21:510–27.
   doi:10.1002/aic.690210313.

.. _ThermoNet: http://www.wiley.com/college/moran/CL_0471465704_S/user/
