from sympy import *
import random, time, math
import numpy as np

p = 71 #Field size for ec choose p st mod(p,4)=3
q = p**2-1 #Field size for program r1cs qap
k = 2 #ec y**2=x**3+a*x+b
a = 1
b = 0
inf = float('inf')

Roots_of_Z = [1,2,3,4]
NumGates = 4
NumWires = 11
Ind_pub = [1,11]
Ind_pri = [2,3,4,5,6,7,8,9,10]
x = symbols('x')
Ax = [0,
        5039*x**3 + 9*x**2 + 5014*x + 24,0,x**3 + 5034*x**2 + 11*x + 5034,\
     3*x**3 + 5016*x**2 + 57*x + 5004,0, x**3 + 5034*x**2 + 11*x + 5034,\
     5037*x**3 + 21*x**2 + 4998*x + 24,0,x**3 + 5034*x**2 + 11*x + 5034,0]
Bx = [x**3 + 5034*x**2 + 11*x + 5034,0,5039*x**3 + 9*x**2 + 5014*x + 24,0,0,\
      3*x**3 + 5016*x**2 + 57*x + 5004,0,0,5037*x**3 + 21*x**2 + 4998*x + 24,0,0]
Cx = [0,0,0, 5034*x**3 + 54*x**2 + 4884*x + 144,0,0,18*x**3 + 4896*x**2 + 342*x + 4824,\
     0,0,5022*x**3 + 126*x**2 + 4788*x + 144,6*x**3 + 5004*x**2 + 66*x + 5004]
Zx = (x-1)*(x-2)*(x-3)*(x-4)
Hx = 5002*x**2 + 19*x + 4971
R = [1,1,3,3,2,4,8,3,4,12,23]
Ap=np.array([[0,0,0,0],
       [24,5014,9,5039],
       [0,0,0,0],
       [5034,11,5034,1],
       [5004,57,5016,3],
       [0,0,0,0],
       [5034,11,5034,1],
       [24,4998,21,5037],
       [0,0,0,0],
       [5034,11,5034,1],
       [0,0,0,0]])
Bp=np.array([[5034,11,5034,1],
       [0,0,0,0],
       [24,5014,9,5039],
       [0,0,0,0],
       [0,0,0,0],
       [5004,57,5016,3],
       [0,0,0,0],
       [0,0,0,0],
       [24,4998,21,5037],
       [0,0,0,0],
       [0,0,0,0]])
Hp = [4971,19,5002]
Px = [5002,399,3451,3255,763,3906,3384]

print('Program to prove: \n x1*x4 + x2*x5 + x3*x6 \n')
print('With input: \n 1 3 2 4 3 4 \n')
print('The number of wires: {wire} \n'.format(wire=NumWires))
print('The number of multiplication gates: {gate} \n'.format(gate=NumGates))
print('Restoration polynomial: ')
print('{a}x^6+{b}x^5+{c}x^4+{d}x^3+{e}x^2+{f}x+{g}'.format(a=Px[0],b=Px[1],\
c=Px[2],d=Px[3],e=Px[4],f=Px[5],g=Px[6]))
print('Quotient polynoomial: ')
print('{a}x^2+{b}x+{c}'.format(a=Hp[0],b=Hp[1],c=Hp[2]))
print('Target polynomial: ')
print(Zx)
print(' ')



#Elliptic curve generation
Curve = {'q':p, 'k':1, 'a':a, 'b':b, 'x':0, 'y':0, 'order':0, 'gen':0}
Curve_gen = Curve.copy()
Curve['gen'] = Curve_gen

Group1 = EC_points(Curve)
Group2 = EC_extend(Group1)
Orders = [grp['order'] for grp in Group1]

#n is order of EC curve group
n = max(Orders)
ind = Orders.index(n)
#g is a generator of EC curve group
g = Group1[ind].copy()
#h is a generator of twisted curve group
h = Group2[ind].copy()
del g['gen']
del h['gen']
for i in range(len(Group1)):
    Group1[i]['gen'] = g.copy()
    Group2[i]['gen'] = h.copy()
PointInf1 = g.copy()
PointInf1['x'] = inf
PointInf1['y'] = inf
PointInf2 = h.copy()
PointInf2['x'] = inf
PointInf2['y'] = inf



#Groth16 Setup
start = time.time()
alpha = 1+random.randrange(1,q-1)
beta = 1+random.randrange(1,q-1)
gamma = 1+random.randrange(1,q-1)
while math.gcd(gamma,q)!=1:
    gamma = 1+random.randrange(1,q-1)
delta = 1+random.randrange(1,q-1)
while math.gcd(delta,q)!=1:
    delta = 1+random.randrange(1,q-1)
x_val = 1
while x_val in Roots_of_Z:
    x_val = random.randrange(q)
#tau is secret
tau = [alpha, beta, gamma, delta, x_val]

Ax_val = []
for i in range(len(Ax)):
    if type(Ax[i])==int:
        Ax_temp = Ax[i]
    else:
        Ax_temp = int(Ax[i].subs(x,x_val))%q
    Ax_val.append(Ax_temp)
Bx_val = []
for i in range(len(Bx)):
    if type(Bx[i])==int:
        Bx_temp = Bx[i]
    else:
        Bx_temp = int(Bx[i].subs(x,x_val))%q
    Bx_val.append(Bx_temp)
Cx_val = []
for i in range(len(Cx)):
    if type(Cx[i])==int:
        Cx_temp = Cx[i]
    else:
        Cx_temp = int(Cx[i].subs(x,x_val))%q
    Cx_val.append(Cx_temp)
Zx_val = Zx.subs(x,x_val)%q
Hx_val = Hx.subs(x,x_val)%q

sigma1_1 = [EC_pmult(alpha,g),EC_pmult(beta,g),EC_pmult(delta,g)]

sigma1_2 = []
for i in range(NumGates):
    sigma1_2.append(0)
    if i==NumGates-1:
        sigma1_2[i] = g
for i in range(NumGates):
    val = (x_val**i)%q
    sigma1_2[i] = EC_pmult(val,g)

sigma1_3 = []
for i in range(NumWires):
    sigma1_3.append(0)
    if i==NumWires-1:
        sigma1_3[i] = g
VAL = [0 for i in range(NumWires)]
for i in range(len(Ind_pub)):
    VAL[Ind_pub[i]-1] = ((beta*Ax_val[Ind_pub[i]-1]\
    +alpha*Bx_val[Ind_pub[i]-1]+Cx_val[Ind_pub[i]-1])*MODinv(gamma,q).real)%q
    sigma1_3[i] = EC_pmult(VAL[Ind_pub[i]-1],g)

sigma1_4 = []
for i in range(NumWires):
    sigma1_4.append(0)
    if i==NumWires-1:
        sigma1_4[i] = g
for i in range(len(Ind_pri)):
    val = ((beta*Ax_val[Ind_pri[i]-1]+alpha*Bx_val[Ind_pri[i]-1]\
    +Cx_val[Ind_pri[i]-1])*MODinv(delta,q).real)%q
    sigma1_4[i] = EC_pmult(val,g)

sigma1_5 = []
for i in range(NumGates-1):
    sigma1_5.append(0)
    if i==NumGates-2:
        sigma1_5[i] = g
for i in range(NumGates-1):
    val = ((x_val**i)*MODinv(delta,q)*Zx_val)%q
    sigma1_5[i] = EC_pmult(val,g)

sigma2_1 = [EC_pmult(beta,h),EC_pmult(gamma,h),EC_pmult(delta,h)]
sigma2_2 = []
for i in range(NumGates):
    sigma2_2.append(0)
    if i==NumGates-1:
        sigma2_2[i] = h
for i in range(NumGates):
    val = (x_val**i)%q
    sigma2_2[i] = EC_pmult(val,h)


print('Groth16 Setup result:')
time.time()-start
print('CRS has been made')
Length_CRS = len(sigma1_1)+len(sigma1_2)+len(sigma1_3)\
+len(sigma1_4)+len(sigma1_5)+len(sigma2_1)+len(sigma2_2)
print('Length of CRS is:'+str(Length_CRS)+'(elliptic curve points)')

#Debugging code
RAx = 0
RBx = 0
RCx = 0
for i in range(len(R)):
    RAx += R[i]*Ax_val[i]
    RBx += R[i]*Bx_val[i]
    RCx += R[i]*Cx_val[i]
if (RAx*RBx-RCx)%q==(Zx_val*Hx_val)%q:
    print('CRS is valid')
else:
    print('CRS is invalid')
print(' ')


#Groth16 Prove
start = time.time()
#Generate r and s which are secret
r = random.randrange(q)
s = random.randrange(q)
#Make proof
Proof_A = sigma1_1[0].copy()
for i in range(NumWires):
    temp = PointInf1.copy()
    for j in range(NumGates):
        temp = EC_add(temp,EC_pmult(Ap[i,j],sigma1_2[j]))
    Proof_A = EC_add(Proof_A,EC_pmult(R[i],temp))
Proof_A = EC_add(Proof_A,EC_pmult(r,sigma1_1[2]))

Proof_B = sigma2_1[0].copy()
for i in range(NumWires):
    temp = PointInf2.copy()
    for j in range(NumGates):
        temp = EC_add(temp,EC_pmult(Bp[i,j],sigma2_2[j]))
    Proof_B = EC_add(Proof_B,EC_pmult(R[i],temp))
Proof_B = EC_add(Proof_B,EC_pmult(s,sigma2_1[2]))

temp_Proof_B = sigma1_1[1].copy()
for i in range(NumWires):
    temp = PointInf1.copy()
    for j in range(NumGates):
        temp = EC_add(temp,EC_pmult(Bp[i,j],sigma1_2[j]))
    temp_Proof_B = EC_add(temp_Proof_B,EC_pmult(R[i],temp))
temp_Proof_B = EC_add(temp_Proof_B,EC_pmult(s,sigma1_1[2]))

Proof_C=EC_add(EC_pmult(s,Proof_A),EC_pmult(r,temp_Proof_B))
Proof_C=EC_add(Proof_C,EC_inv(EC_pmult(r,EC_pmult(s,sigma1_1[2]))))
for i in range(len(Ind_pri)):
    Proof_C=EC_add(Proof_C,EC_pmult(R[Ind_pri[i]-1],sigma1_4[i]))
for i in range(NumGates-1):
    Proof_C=EC_add(Proof_C,EC_pmult(Hp[i],sigma1_5[i]))
proof = [Proof_A,Proof_B,Proof_C]
print('Groth16 Prove result:')
time.time()-start
print('Proof has been made')
print('Length of proof is:'+str(len(proof))+'(elliptic curve points)')

#Debuging code
A = (alpha+RAx+r*delta)%q
B = (beta+RBx+s*delta)%q
temp_C=0
for i in Ind_pri:
    temp_C+=R[i-1]*(beta*Ax_val[i-1]+alpha*Bx_val[i-1]+Cx_val[i-1])
C = (MODinv(delta,q)*(temp_C+Hx_val*Zx_val)\
+A*s+B*r+((-r*s*delta)%q))%q

lhs = (A*B)%q
temp_rhs=0
for i in Ind_pub:
    temp_rhs+=R[i-1]*VAL[i-1]
rhs = (alpha*beta+gamma*temp_rhs+C*delta)%q

proofcheckflag = 1
proofcheckflag = proofcheckflag*Peq(Proof_A,EC_pmult(A,g))\
*Peq(Proof_B,EC_pmult(B,h))*Peq(Proof_C,EC_pmult(C,g))
if proofcheckflag==1 and lhs==rhs:
    print('Proof is complete')
else:
    print('Proof is incomplete')
print(' ')
print('Proof')
print('{a}\n{b}\n{c}\n'.format(a=proof[0],b=proof[1],c=proof[2]))

#Groth16 Verify
start = time.time()
LHS = cmod(weil(proof[0],proof[1],Group1,Group2),p)
RHS = 1
RHS = cmod(RHS*weil(sigma1_1[0],sigma2_1[0],Group1,Group2),p)
temp = PointInf1.copy()

for i in range(len(Ind_pub)):
    temp = EC_add(temp,EC_pmult(R[Ind_pub[i]-1],sigma1_3[i]))
RHS = cmod(RHS*weil(temp,sigma2_1[1],Group1,Group2),p)
RHS = cmod(RHS*weil(proof[2],sigma2_1[2],Group1,Group2),p)

VfyResult = LHS==RHS
print('Groth16 Verify result:')
time.time()-start
if VfyResult==1:
    print('Verification Success')
else:
    print('Verification Failure')
