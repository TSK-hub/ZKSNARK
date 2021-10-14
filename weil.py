import random

def weil(Point1,Point2,G1,G2):
    a = Point1['a']
    p = Point1['q']
    P = Point1.copy()
    Q = Point2.copy()
    n = G1[1]['gen']['order']
    inf = float('inf')
    O = [inf,inf]

    bstring = format(n,'b')
    m = len(bstring)
    RP = P.copy()
    PpRP = EC_add(P,RP)
    while Peq(RP,PpRP) or (PpRP['x']==inf and PpRP['y']==inf):
        idx = random.randrange(n)
        RP = G1[idx-1].copy()
        PpRP = EC_add(P,RP).copy()
        RP = P.copy()
        PpRP=EC_add(P,RP)
    RQ = Q.copy()
    QpRQ = EC_add(Q,RQ)
    while Peq(RQ,QpRQ) or (QpRQ['x']==inf and QpRQ['y']==inf):
        idx = random.randrange(n)
        RQ = G2[idx-1].copy()
        QpRQ = EC_add(Q,RQ).copy()
        RQ=Q.copy()
        QpRQ = EC_add(Q,RQ)

    numer = cmod(Miller(n,Q,P,P),p)*cmod(Miller(n,P,Q,Q),p)
    denom=cmod(Miller(n,Q,P,P),p)*cmod(Miller(n,P,Q,Q),p)
    Weil=cmod(cmod(numer,p)*MODinv(denom,p),p)
    return Weil
