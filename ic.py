import numpy as np

# Set intial condition for burger's equation

def dflu1(x):
    return 0.0

def dflu2(x):
    if x < 0.0:
        return 0.9
    else:
        return 0.2
dflu2 = np.vectorize(dflu2)

def smooth(x):
    return np.sin(np.pi*x)
    
def composite(x):
    f = np.empty_like(x)
    for i,y in enumerate(x):
        if y > -0.8 and y < -0.6:
            f[i] = np.exp(-np.log(2.0)*(y+0.7)**2/0.0009);
        elif y > -0.4 and y < -0.2:
            f[i] = 1.0;
        elif y > 0.0 and y < 0.2:
            f[i] = 1.0 - abs(10.0*(y-0.1));
        elif y > 0.4 and y < 0.6:
            f[i] = np.sqrt(1.0 - 100.0*(y-0.5)**2);
        else:
            f[i] = 0.0;
    return f
# Chiarello 2020 paer
def smooth_c(x):  
    return 0.5+(0.4*np.sin(np.pi*x))

def shock(x):
    
    if x < 0.0:
        u = 1.0
    else:
        u = 0.0
    return u
shock = np.vectorize(shock)

def hat(x):
    if (x> -0.75 and x <-0.25):
        u = 1.0
    else:
        u = 0.0
    return u
hat = np.vectorize(hat)

def gauss(x):
    return 0.8 + np.exp(-((x - 0.5)**2) / (0.05**2))
gauss = np.vectorize(gauss)
def triangle(x):
    if 0.0 <= x < 0.3:
        return 0.0  # linear ramp up from 0 to 1
    elif 0.3 <= x and x <= 0.5:
        return (x-0.3)/0.2
    elif 0.5 < x and x <= 0.7:
        return (0.7 - x) / 0.2  # linear ramp down from 1 to 0
    else:
        return 0.0
triangle = np.vectorize(triangle)
    

def buckley1(x):
    if x>= (-1.0/2.0) and x <= 0:
        return 1.0
    else:
        return 0
buckley1 = np.vectorize(buckley1)

# Rarefaction without sonic point
def rare1(x):
    if x < 0.:
        u = -1
    else:
        u = 1
    return u
rare1 = np.vectorize(rare1)

# Rarefaction with sonic point
def rare(x):
    #u = np.zeros(len(x))
    #for i in range(len(x)):
    if x < 0.5:
        u = -0.5
    else:
        u = 1.0
    return u
rare = np.vectorize(rare)

def rare3(x):
    #u = np.zeros(len(x))
    #for i in range(len(x)):
    if x < 0.0:
        u = 0.0
    elif x >=0.0 and x<=1.0:
        u = 1.0
    else:
        u = 0.0    
    return u
rare3 = np.vectorize(rare3)

def expo(x):
    return 1.0 + np.exp(-100*(x-0.25)**2)
