import math,cmath

def eval_f1(T,P,R):
    #eval f1 at the point T
    # div(f1)=(P+R)-(R)-(P)+(O)
    inf = float('inf')
    p = T['q']
    O = [inf,inf]
    M = EC_add(P,R).copy()
    g1 = EC_line(P,R).copy()
    if M['x']==inf and M['y']==inf:
        g2 = [0,0,1]
    else:
        g2 = EC_line(M,EC_inv(M)).copy()
    out = cmod(eval_line(g2,T)*MODinv(eval_line(g1,T),p),p)
    return out
