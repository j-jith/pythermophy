from __future__ import division, print_function
import yaml

class Fluid(object):
    """
    Class for representing the fluid you want to study.

    :param name: name of the fluid
    :type name: str
    :param M: molar mass of the fluid (kg/mol)
    :type M: float
    :param Tc: critical temperature of the fluid (K)
    :type Tc: float
    :param pc: critical pressure of the fluid (Pa)
    :type pc: float
    :param acentric: acentric factor of the fluid (dimensionless)
    :type acentric: float
    :param cp_coeffs: coefficients of the ideal gas isobaric heat capacity
        polynomial [a,b,c,d] :math:`c_p = a + bT + cT^2 + dT^3` (:math:`c_p` in
        J/mol/K)
    :type cp_coeffs: list(float)
    """

    def __init__(self, name, M, Tc, pc, acentric, cp_coeffs):
        """
        Creates a new :class:`Fluid` instance.
        """

        self.name = name
        self.molar_mass = M
        self.T_crit = Tc
        self.p_crit = pc
        self.ideal_cp_coeffs = cp_coeffs
        self.acentric = acentric

    def is_valid(self):
        """
        Checks if the :class:`Fluid` instance is valid by ensuring that all the
        parameters of the object are of the correct type.

        :return: True if valid. False otherwise.
        :rtype: bool
        """

        invalid = []

        if not isinstance(self.name, basestring):
            invalid.append('name should be a string')
        if not isinstance(self.molar_mass, int) and not isinstance(self.molar_mass, float):
            invalid.append('molar_mass should be int/float')
        if not isinstance(self.T_crit, int) and not isinstance(self.T_crit, float):
            invalid.append('T_crit should be int/float')
        if not isinstance(self.p_crit, int) and not isinstance(self.p_crit, float):
            invalid.append('p_crit should be int/float')
        if not isinstance(self.acentric, int) and not isinstance(self.acentric, float):
            invalid.append('acentric should be int/float')

        for item in self.ideal_cp_coeffs:
            if not isinstance(item, int) and not isinstance(item, float):
                invalid.append('ideal_cp_coeffs should be a list of int/float')
                break

        if invalid:
            return (False, invalid)
        else:
            return (True, None)


    @classmethod
    def init_from_file(cls, filename):
        """
        Creates a new :class:`Fluid` instance by reading parameters from an
        input file.

        The input file must be a ``yaml`` file in the following format:

        .. code:: yaml

            name: CO2   # Name of the fluid
            molar_mass: 44.01e-3    # Molar mass [kg/mol]
            T_crit: 304.25  # Critical temperature [K]
            p_crit: 7.38e+6 # Critical pressure [Pa]
            acentric: 0.228 # Acentric factor [dimensionless]
            ideal_cp_coeffs: [22.26, 5.981e-2, -3.501e-5, 7.469e-9] # Cp = a + b*T + c*T^2 + d*T^3 [J/mol/K]

        :param filename: input file
        :type filename: str
        """

        try:
            with open(filename, 'r') as fileh:
                data = yaml.load(fileh)
                #print(data)
                try:
                    name = data['name']
                    M = data['molar_mass']
                    Tc = data['T_crit']
                    pc = data['p_crit']
                    acentric = data['acentric']
                    cp_coeffs = data['ideal_cp_coeffs']

                    inst = cls(name, M, Tc, pc, acentric, cp_coeffs)
                    valid = inst.is_valid()
                    if not valid[0]:
                        raise TypeError('These items are not of the right type in the file: ' + ', '.join(valid[1]))
                    else:
                        return inst

                except KeyError as e:
                    print('Error: File does not contain data in the correct format.', e)
                except TypeError as e:
                    print('Error:', e)
        except IOError as e:
            print('Cannot read file.', e)


    def save_to_file(self, filename):
        """
        Saves the :class:`Fluid` instance to a yaml file which can later be fed to
        :py:func:`init_from_file`.

        :param filename: output file
        :type filename: str
        """

        with open(filename, 'w') as fileh:
            yaml.dump(self.__dict__, fileh)
