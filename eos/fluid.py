from __future__ import division, print_function
import yaml

class Fluid(object):

    def __init__(self, name, M, Tc, pc, acentric, cp_coeffs):

        self.name = name
        self.molar_mass = M
        self.T_crit = Tc
        self.p_crit = pc
        self.ideal_cp_coeffs = cp_coeffs
        self.acentric = acentric

    def is_valid(self):

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

        with open(filename, 'w') as fileh:
            yaml.dump(self.__dict__, fileh)
