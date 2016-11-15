from sympy import *

init_printing()

p, T, R, d, D = symbols('p,T,R,d,D')
V = Function('V')(T)
a = Function('a')(T)
Z = Function('Z')(T)

# Cp
h = Function('h')(T)

h1 = (2*V + d + D)/(2*V +d - D)

H = p*V - R*T + (T*diff(a, T) - a)/D * log(h)
H = H.subs({V:Z*R*T/p})

cp = diff(H, T)

h1 = h1.subs({V:Z*R*T/p})
dh_T = diff(h1, T).simplify()

# Cv
cv = diff((T*diff(a, T) - a)/D, T) * log(h)

