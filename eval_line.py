def eval_line(g,Point):
    #eval point P=(x,y)
    p = Point['q']
    P = [Point['x'],Point['y']]
    xp = P[0]
    yp = P[1]
    return cmod(g[0]*yp+g[1]*xp+g[2],p)
