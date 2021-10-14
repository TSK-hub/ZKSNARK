import math,cmath

def Miller(n, T, P, R):
    inf = float('inf')
    Z = P.copy()
    Z['x'] = inf
    Z['y'] = inf
    V = 1
    k = 0
    f1 = eval_f1(T, P, R)
    bstring = format(n, 'b')
    m = len(bstring)

    for i in range(m):
        if i < m-1:
            V = Algo_D(T, V, V, Z, Z)
            Z = EC_pmult(2, Z).copy()
            k *=2
        if int(bstring[m-i-1]) == 1:
            V = Algo_D(T, V, f1, Z, P)
            Z = EC_add(Z, P).copy()
            k += 1
        

    return V
