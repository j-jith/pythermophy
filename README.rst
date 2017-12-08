**pythermophy** is a python module for predicting various thermophysical
properties of real fluids. The fluids can be currently  modelled using the
following equations of state:

- Ideal Gas (IG)
- Redlich-Kwong (RK)
- Soave-Redlich-Kwong (SRK)
- Peng-Robinson (PR)
- Lee-Kesler (LK)

At the moment, the following thermophysical properties can be predicted as a
function of temperature ``T`` (K) and pressure ``p`` (Pa):

- Density (kg/m^3)
- Compressibility factor
- Speed of sound (m/s)
- Isobaric specific heat capacity (real and ideal) (J/mol/K)
- Isochoric specific heat capacity (real and ideal) (J/mol/K)
- Departure function of specific heat capacities (J/mol/K)
- Adiabatic index
- Isothermal compressibility (1/Pa)

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

Documentation
=============

For detailed documentation, please visit https://j-jith.github.io/pythermophy/
