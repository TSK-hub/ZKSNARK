import cmath

def EC_inv(Point):
    inf = float('inf')
    p = Point['q']
    P = [Point['x'],Point['y']]
    Point2 = Point.copy()
    O = [inf,inf]
    if P==O:
        Q = P.copy()
    elif P[1]==0:
        if Point['k']>1:
            Q = [P[0],p*cmath.sqrt(-1)]
        else:
            Q = [P[0],p]
    else:
        Q = [P[0], cmod(-P[1],p)]
    Point2['x'] = Q[0]
    Point2['y'] = Q[1]
    return Point2
