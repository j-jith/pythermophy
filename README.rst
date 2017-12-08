.. default-role:: math

**pythermophy** is a python module for predicting various thermophysical
properties of real fluids. The fluids can be currently  modelled using the
following equations of state:

- Ideal Gas (IG)
- Redlich-Kwong (RK)
- Soave-Redlich-Kwong (SRK)
- Peng-Robinson (PR)
- Lee-Kesler (LK)

At the moment, the following thermophysical properties can be predicted as a
function of temperature `T` (K) and pressure `p` (Pa):

- Density, `\rho` (kg/m^3)
- Compressibility factor, `Z`
- Speed of sound, `c` (m/s)
- Isobaric specific heat capacity (real and ideal), `c_p` (J/mol/K)
- Isochoric specific heat capacity (real and ideal), `c_v` (J/mol/K)
- Departure function of specific heat capacities (J/mol/K)
- Adiabatic index, `\gamma`
- Isothermal compressibility, `\chi_T` (1/Pa)

Installation
============

To install pythermophy, it is recommended to clone the git repository and use
the provided ``setup.py`` as follows.

.. code::

    $ git clone https://github.com/j-jith/pythermophy
    $ cd pythermophy
    $ python setup.py install

The module was written and tested on Python 3, and therefore it is recommended
to use pythermophy with Python 3.

