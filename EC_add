import cmath

def EC_add(Point1,Point2):
    if Point1['k']!=Point2['k']:
        raise ValueError('mismatched degrees')
    R = Point1.copy()
    a = Point1['a']
    p = Point1['q']
    k = Point1['k']

    P = [Point1['x'],Point1['y']]
    Q = [Point2['x'],Point2['y']]
    nan = float('nan')
    inf = float('inf')
    O = [inf,inf]

    if P==O and Q==O:
        R['x'] = inf
        R['y'] = inf
        return R
    else:
        if P==O:
            R = Point2.copy()
            return R
        elif Q==O:
            R = Point1.copy()
            return R
        
    if cmath.isnan(P[0]*Q[0]) and cmath.isnan(P[1]*Q[1]):
        R['x'] = nan
        R['y'] = nan
        return R
    

    if k>1:
        P = [cmod(-(Point1['x']),p), Point1['y']/cmath.sqrt(-1)]
        Q = [cmod(-(Point2['x']),p), Point2['y']/cmath.sqrt(-1)]
    else:
        P = [Point1['x'], Point1['y']]
        Q = [Point2['x'], Point2['y']]

    xp = P[0]
    xq = Q[0]
    yp = P[1]
    yq = Q[1]
    if xp!=xq:
        s = cmod((yq-yp)*MODinv((xq-xp),p),p)
    else:
        if yp!=yq:
            R['x'] = inf
            R['y'] = inf
            return R
        else:
            if yp==0:
                R['x'] = inf
                R['y'] = inf
                return R
            else:
                s = cmod((3*(xp**2)+a)*MODinv((2*yp),p),p)
    xr = cmod(s**2-xp-xq,p)
    yr = cmod(s*(xp-xr)-yp,p)

    if k>1:
        R['x'] = cmod(-xr,p)
        R['y'] = yr*cmath.sqrt(-1)
        return R
    else:
        R['x'] = xr
        R['y'] = yr
        return R
