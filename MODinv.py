import cmath, math

def MODinv(a,n):
    #the inverse of a modulo n
    inf = float('inf')
    nan = float('nan')
    if a==0:
        raise ValueError('divide by zero')
    if a.imag==0:
        mult = 1
        b = a
    else:
        mult = a.real-a.imag*cmath.sqrt(-1)
        b = a.real**2+a.imag**2
    if b==inf or b==-inf:
        raise ValueError('divide by inf or -inf')
    while b.real<0:
        b+=n
    t = 0
    r = n
    newt = 1
    newr = b.real
    while newr!=0:
        Quotient = math.floor(r/newr)
        oldt = t
        t = newt
        newt = oldt-Quotient*newt
        oldr = r
        r = newr
        newr = oldr-Quotient*newr
        if cmath.isnan(newr):
            raise ValueError('No inverse exists')
    if r>1:
        aINV=nan
        return aINV
    if t<0:
        t+=n
    aINV = cmod(t*mult,n)
    return aINV
