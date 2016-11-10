from sympy import *

init_printing()

p, T = symbols('p, T')

A = Function('A')(p, T)
B = Function('B')(p, T)
C = Function('C')(p, T)
Z = Function('Z')(p, T)

eq = Z**3 + A*Z**2 + B*Z + C

deq_T = diff(eq, T)
dZ_T = solve(deq_T, diff(Z, T))[0]

deq_p = diff(eq, p)
dZ_p = solve(deq_p, diff(Z, p))[0]
