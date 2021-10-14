def EC_points(Curve):
    p = Curve['q']
    a = Curve['a']
    b = Curve['b']
    x = [i for i in range(p)]
    y = [i for i in range(p)]
    idxs = []
    for i in range(len(x)):
        for j in range(len(y)):
            z = (j**2)%p - (i**3+a*i+b)%p
            if z==0:
                idxs.append([i,j])
    size = len(idxs)
    temp = Curve.copy()
    Group = []
    for i in range(size):
        ref = idxs[i]
        temp['k'] = 1
        temp['x'] = ref[0]
        temp['y'] = ref[1]
        Group.append(EC_order(temp))

    return Group
