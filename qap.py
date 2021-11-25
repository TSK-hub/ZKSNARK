import numpy as np
import sympy as sp
import math

p=71
q=p**2-1

A=np.array([[0,0,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,0,0],\
          [0,0,0,1],[0,0,1,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]])
B=np.array([[0,0,0,1],[0,0,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,1,0,0],\
          [0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0],[0,0,0,0]])
C=np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0],\
          [0,1,0,0],[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]])
R=[1,3,20,60,4,30,120,3,35,105,285]

x=sp.Symbol('x')

NumGates=len(A[0])
NumWires=len(A)

Roots_of_Z=[i+1 for i in range(NumGates)]
Zx=1
for i in range(NumGates):
    Zx*=(x-Roots_of_Z[i])
Zp=[int(sp.expand(Zx).coeff(x,i)) for i in range(NumGates+1)]

Ap=np.zeros((NumWires,NumGates))
Bp=Ap.copy()
Cp=Bp.copy()

VANDER=np.flip(np.vander(Roots_of_Z),1)
for i in range(NumWires):
    y_vec=A[i].copy()
    polycoeffs=np.linalg.solve(VANDER,y_vec)
    Ap[i]=polycoeffs.copy()
    
for i in range(NumWires):
    y_vec=B[i].copy()
    polycoeffs=np.linalg.solve(VANDER,y_vec)
    Bp[i]=polycoeffs.copy()
    
for i in range(NumWires):
    y_vec=C[i].copy()
    polycoeffs=np.linalg.solve(VANDER,y_vec)
    Cp[i]=polycoeffs.copy()
    
denoms=[]
for i in range(NumWires):
    for j in range(NumGates):
        denom=int(sp.fraction(sp.nsimplify(Bp[i][j]))[1])
        if denom not in denoms:
            denoms.append(denom)
for i in range(NumWires):
    for j in range(NumGates):
        denom=int(sp.fraction(sp.nsimplify(Bp[i][j]))[1])
        if denom not in denoms:
            denoms.append(denom)
for i in range(NumWires):
    for j in range(NumGates):
        denom=int(sp.fraction(sp.nsimplify(Cp[i][j]))[1])
        if denom not in denoms:
            denoms.append(denom)
denoms.sort()
mul=1
for i in range(len(denoms)):
    GCD=math.gcd(int(mul),denoms[i])
    mul*=denoms[i]
    mul/=GCD
    
Ap=np.round(Ap*mul)
Bp=np.round(Bp*mul)
Cp=np.round(Cp*mul*mul)

Ap=rmod(Ap,q)
Bp=rmod(Bp,q)
Cp=rmod(Cp,q)

Zp=rmod([Zp],q)
Zp=Zp[0]
R=rmod([R],q)
R=R[0]

xv=[]
for k in range(NumGates):
    xv.append(x**k)

Ax=Ap.dot(xv)
Bx=Bp.dot(xv)
Cx=Cp.dot(xv)

Px=Ax.dot(R)*Bx.dot(R)-Cx.dot(R)
Px=sp.simplify(Px)
Px=pmod(Px,q)

checksum=0
for i in Roots_of_Z:
    checksum+=float(Px.subs(x,i))%q
if checksum!=0:
    print('invalid P(x), check the interpolation')

Hx=390*x**2+189*x+4005
Hp=[4005,189,390]
Ax=list(Ax)
Bx=list(Bx)
Cx=list(Cx)
R=list(R)
