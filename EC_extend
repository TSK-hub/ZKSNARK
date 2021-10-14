import cmath

def EC_extend(Group1):
    #work only for curve of the form y^2=x^3+ax w/ p=3mod4
    #assumed k=2
    p = Group1[0]['q']
    Group2 = []
    size = len(Group1)
    for i in range(size):
        Group2.append(Group1[i].copy())
        Group2[i]['k'] = 2
        Group2[i]['x'] = cmod(((-1)*Group1[i]['x']),p)
        Group2[i]['y'] = cmath.sqrt(-1)*Group1[i]['y']
    return Group2
