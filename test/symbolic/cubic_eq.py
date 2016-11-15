from sympy import *

init_printing()

p, V, T, R, b, c, d, e = symbols('p,V,T,R,b,c,d,e')
a = Function('a')(T)
Z = symbols('Z')

# general cubic equation of state
f = p - R*T/(V-b+c) + a/(V**2+d*V+e)

# multiply equation of state by common denominator
g = f*((V-b+c)*(V**2+V*d+e))

# substitute V=ZRT/p
h = g.subs({V:Z*R*T/p})

# simplify and expand
ff = h.simplify().expand()

# multiply throughout to make coeff. of Z unity
gg = ff*p**2/R**3/T**3

# expand and collect the expression in terms of Z
hh = gg.expand().collect(Z)

# create a polynomial from the expression
eq = Poly(gg, Z)

# obtain the coefficients of the polynomial
coeffs = eq.coeffs()
