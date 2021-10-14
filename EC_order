import cmath

def EC_order(Point):
    modified_Point = Point.copy()
    inf = float('inf')
    O = [inf,inf]
    a = Point['a']
    p = Point['q']
    k = Point['k']
    if k>1:
        P = [-(Point['x'])%p, Point['y']/cmath.sqrt(-1)]
    else:
        P = [Point['x'], Point['y']]

    if P[0]==inf or P[1]==inf:
        n = 1
        modified_Point['order'] = n
        return modified_Point

    Q = EC_add(Point,Point)
    if Q['x']==inf or Q['y']==inf:
        n = 2
        modified_Point['order'] = n
        return modified_Point

    for n in range(2,p+1):
        Q = EC_add(Point,Q)
        if Q['x']==inf or Q['y']==inf:
            modified_Point['order'] = n+1
            break
    return modified_Point
