import sympy as sp

def pmod(poly,q):
    newpoly=0
    deg=int(sp.degree(poly))
    for i in range(deg+1):
        coef=poly.coeff(x,i)%q
        newpoly+=coef*(x**i)
    return newpoly
