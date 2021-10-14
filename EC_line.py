import cmath

def EC_line(Point1,Point2):
    #output: g(x,y)=cy+dx+e
    if Point1['k']!=Point2['k']:
        raise ValueError('mismatched degrees')
    a = Point1['a']
    p = Point1['q']
    k = Point1['k']
    if k>1:
        P = [cmod(-Point1['x'],p), Point1['y']/cmath.sqrt(-1)]
        Q = [cmod(-Point2['x'],p), Point2['y']/cmath.sqrt(-1)]
    else:
        P = [Point1['x'],Point1['y']]
        Q = [Point2['x'],Point2['y']]

    xp = P[0]
    xq = Q[0]
    yp = P[1]
    yq = Q[1]
    inf = float('inf')

    O = [inf,inf]
    if P==O and Q==O:
        raise ValueError('Both input points are Point at Infinity')


    if P==Q and yp!=0:
        #tangent
        s = cmod((3*(xp**2)+a)*MODinv(2*yp,p),p)
        c = 1
        d = cmod(-s,p)
        e = cmod(s*xp-yp,p)
    elif xp==xq:
        #vertical
        c = 0
        d = 1
        e = cmod(-xp,p)
    elif P==O and Q!=O:
        #vertical
        c = 0
        d = 1
        e = cmod(-xq,p)
    elif P!=O and Q==O:
        #vertical
        c = 0
        d = 1
        e = cmod(-xp,p)
    else:
        s = cmod((yq-yp)*MODinv((xq-xp),p),p)
        c = 1
        d = cmod(-s,p)
        e = cmod(s*xp-yp,p)
    if k>1:
        c = c/cmath.sqrt(-1)
        d = cmod(-d,p)
    return [c,d,e]
