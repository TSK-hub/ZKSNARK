import cmath

def cmod(val,q):
    inf = float('inf')
    if val.real*val.imag==inf or val.real*val.imag==-inf:
        print('Warn : modular of inf or -inf')
        out = val
        return out
    else:
        out = val.real%q+val.imag%q*cmath.sqrt(-1)
        return out
