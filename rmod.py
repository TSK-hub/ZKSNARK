import numpy as np
import sympy as sp

def rmod(mat,q):
    #Take modular on rational numbers
    Y=np.zeros((len(mat),len(mat[0])))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            x=mat[i][j]
            if round(x)==x:
                y=x%q
            else:
                (n,d)=sp.fraction(sp.nsimplify(x))
                y=(int(n)*MODinv(int(d),q))%q
            Y[i][j]=y
    return Y
