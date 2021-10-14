def EC_pmult(p,Point):
    Q = Point.copy()
    inf = float('inf')
    p = p%(Point['q']+1)
    if p==0:
        Q['x'] = inf
        Q['y'] = inf
        return Q
    elif p==1:
        return Q
    Q = EC_add(Point,Point)
    for i in range(int(p)):
        if i==0 or i==1:
            continue
        else:
            Q = EC_add(Point,Q)
    return Q
