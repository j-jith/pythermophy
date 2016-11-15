from sympy import *

init_printing()

# constants
b1, b2, b3, b4 = symbols('b1, b2, b3, b4')
c1, c2, c3, c4 = symbols('c1, c2, c3, c4')
d1, d2 = symbols('d1, d2')
beta, gamma = symbols('beta, gamma')

# p, V, T
pr, Vr, Tr = symbols('pr, Vr, Tr')

Z = Function('Z')(pr, Tr)

# coefficients
B = Function('B')(Tr)
C = Function('C')(Tr)
D = Function('D')(Tr)

eq = 1 + B/Vr + C/Vr**2 + D/Vr**5 + c4/Tr**3/Vr**2 * (beta + gamma/Vr**2) * exp(-gamma/Vr**2) - pr*Vr/Tr

eq = eq.subs({Vr: Tr*Z/pr})

dZ_Tr = solve(diff(eq, Tr), diff(Z, Tr))
