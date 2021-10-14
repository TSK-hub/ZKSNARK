def Algo_D(T, f_b, f_c, Pb, Pc):
    #evaluate f_{b+c} at the point T
    inf = float('inf')
    p = T['q']
    O = [inf,inf]

    if (Pb['x']==inf and Pb['y']==inf) or (Pc['x']==inf and Pc['y']==inf):
        return cmod(f_b*f_c,p)

    Pd = EC_add(Pb,Pc).copy()
    g1 = EC_line(Pb,Pc).copy()

    if Pd['x']*Pd['y']==0:
        g2 = [0,0,1]
    else:
        g2 = EC_line(Pd,EC_inv(Pd)).copy()

    return cmod(f_b*f_c*eval_line(g1,T)*MODinv(eval_line(g2,T),p),p)
    
