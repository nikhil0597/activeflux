import numpy as np
#set domain
xmin, xmax = 0.0, 1.0
# f = u
def flux(u):
    c =1.0
    return c*u

# works for any initial data    
def uexact(x, t, u0):
    c = 1.0                
    return  u0(x-c*t)
    #return u0(x)

