"""
Solve u_t + cu_x = 0  for c constant 
Semi-discrete active flux scheme
"""
from xml import dom
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *
from scipy import optimize
from scipy.interpolate import interp1d
from scipy.integrate import quad


# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-nc', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.01)
parser.add_argument('-ic',
                    choices=('gauss', 'triangle', 'smooth', 'composite', 'hat'),
                    help='Initial condition', default='gauss')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)

parser.add_argument('-pde', choices=('linear','burger'),help='PDE', default='linear')
parser.add_argument('-compute_error', choices=('no','yes'),
                    help='Compute error norm', default='no')
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution',
                    default=1)
parser.add_argument('-time_scheme', default = 'euler',
                    help = 'Chosen by degree if unspecified',
                    choices = ('euler','ssprk2', 'ssprk3'))
parser.add_argument('-bc', default = 'periodic',
                    help = 'Choose the boundary condition',
                    choices = ('periodic','dirichlet'))
args = parser.parse_args()


# Select PDE
if args.pde == 'burger':
    from burger import *
elif args.pde == 'linear':
    from linear import *

# constants
Tf    = args.Tf
cfl   = args.cfl
uinit = args.ic
nc    = args.nc
time_scheme = args.time_scheme

if args.ic == 'smooth':
    uinit = smooth
elif args.ic == 'triangle':
    uinit = triangle
elif args.ic == 'hat':
    uinit = hat
elif args.ic == 'gauss':
    uinit = gauss  
elif args.ic == 'composite':
    uinit = composite       
x   = np.zeros(nc)
x_int = np.zeros(nc+1)
h = (xmax - xmin)/nc
  
for i in range(nc):
    x[i] = xmin + i*h + 0.5*h      #cell-center points, nc in number

for i in range(nc+1):
    x_int[i] = xmin + i*h          #interface points, nc+1 in number 

u = np.zeros(nc+2)                 #cell-average solution variable, includes ghost cells
for i in range(nc):                #initialize cell-average solution variable 
    a = x_int[i]
    b = x_int[i+1]
    u[i+1], _ = quad(uinit, a, b)  # integrate initial datum over the cells
    u[i+1] /= h                    # divide by h to get the average

u_int = np.zeros(nc+3)             #interface solution variable, includes ghost cells   
u_int[1:nc+2] = uinit(x_int)       #initialize interface solution variable
res = np.zeros(nc+2)               #residue corresponding to average values
w = np.zeros(nc+3)                 #residue corresponding to interface values


t = 0.0 # time
# plot initial condition
if args.plot_freq >0:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #line1,line2 = ax.plot(x, u[1:nc+1], 'ro',x, u[1:nc+1], 'b')
    line1,line2 = ax.plot(x_int, u_int[1:nc+2], 'ro',x_int, u_int[1:nc+2], 'b')
    #line1 = ax.plot(x, u[2:nc+2], 'ro')
    #line1, = ax.plot(x, u, 'o')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.title('nc='+str(nc)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
    plt.legend(('Approx. solution','Exact solution'))
    #plt.legend(('Solution'))
    #plt.ylim(0.3,0.7)
    plt.grid(True); plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

#Error computation for pointwise values
def compute_error(u1,x1,t):
    error_norm1 = 0.0; error_norm2 = 0.0
    ue = uexact(x1, t , uinit)
    dom_len = xmax - xmin
    error_norm1 = h*np.sum(np.abs(u1-ue))
    error_norm2 = np.sqrt(h*np.sum((u1-ue)**2))
    return error_norm1, error_norm2

#Error computation for average values
def compute_error_avg(u1,x1,t):
    error_norm1 = 0.0; error_norm2 = 0.0
    def ex(x): 
        return uexact(x, t, uinit)    #exact solution at time t
    z = np.zeros(nc)                  #variable for saving cell-averages of exact solution
    for i in range(nc):
        a = x_int[i]
        b = x_int[i+1]
        z[i], _ = quad(ex, a, b)      #integrate interpolated exact solution over the cells
        z[i] /= h                     #divide by h to get the average
    dom_len = xmax - xmin
    error_norm1 = h*np.sum(np.abs(u1-z))
    error_norm2 = np.sqrt(h*np.sum((u1-z)**2))
    return error_norm1, error_norm2




def update_ghost(u1, v1):
    if args.bc == 'periodic':
        # left ghost cell values
        u1[0] = u1[nc]
        v1[0] = v1[nc]    

        # right ghost values
        u1[nc+1] = u1[1]
        v1[nc+2] = v1[2]
    else:
        print('unknown boundary condition')
        exit()
    
# First order Euler forward step
def apply_euler(t,lam, u_old, v_old, u, v, ures, w):
    ts  = t
    update_ghost(u, v)
    ures, w = compute_residual(ts, lam, u, v, ures, w)
    u = u - lam * ures
    v = v - lam*w
    return u, v

def apply_ssprk2(t,lam, u_old, v_old, u, v, ures, w):
    #first stage
    ts  = t
    update_ghost(u,v)
    ures, w = compute_residual(ts, lam, u, v, ures, w)
    u = u - lam * ures
    v = v - lam*w
    #second stage
    ts = t + dt
    update_ghost(u,v)
    ures, w = compute_residual(ts, lam, u, v, ures, w)
    u = 0.5 * u_old + 0.5 *(u - lam * ures)
    v = 0.5 * v_old + 0.5 *(v - lam * w)
    return u, v

def apply_ssprk3(t,lam, u_old, v_old, u, v, ures, w):
    # First stage
    ts = t
    update_ghost(u,v)
    ures, w = compute_residual(ts, lam, u, v, ures, w)
    u = u - lam * ures
    v = v- lam*w

    # Second stage
    ts = t + dt
    update_ghost(u, v)
    ures, w = compute_residual(ts, lam, u, v, ures, w)
    u = (0.75 * u_old) + 0.25*(u - lam * ures)
    v = (0.75 * v_old) + 0.25*(v - lam * w)

    # Third stage
    ts = t + 0.5 * dt  
    update_ghost(u,v)
    ures, w = compute_residual(ts, lam, u, v, ures, w)
    u = (1.0 / 3.0) * u_old + (2.0 / 3.0) * (u - lam * ures)
    v = (1.0 / 3.0) * v_old + (2.0 / 3.0) * (v - lam * w)
    return u, v


def compute_residual(ts, lam, u, v, res, w):
    res[:] = 0.0
    w[:] = 0.0    
    for i in range(1,nc+2):      # face between i-1 and i
        xf = xmin+(i-1)*h        # location of the face
        vl= v[i]
        fl= flux(vl)
        res[i] -= fl
        res[i-1] += fl
        w[i] = (5.0/6.0)*v[i-1] - 3.0*u[i-1] + (4.0/3.0)*v[i] + u[i] - (1.0/6.0)*v[i+1]
    res[0], res[nc+1] = 0.0, 0.0        
    return res, w
time_schemes = {'euler': apply_euler, 'ssprk2' : apply_ssprk2, 'ssprk3' : apply_ssprk3}
t, it = 0.0, 0
while t < Tf:
    if args.pde == 'linear':
        dt= cfl * h
    else:
        print('dt is not set')
        exit()
    lam = dt/h
    if t+dt > Tf:
        dt = Tf - t
        lam = dt/h
    u_old = u    
    v_old = u_int
    u, u_int = time_schemes[time_scheme](t, lam, u_old, v_old, u, u_int, res, w)  # update solution
    t += dt; it += 1 # update time step
    if args.plot_freq >0:
        #ue = uexact(x, t , uinit)
        ue = uexact(x_int, t , uinit)
        #line1.set_ydata(u[1:nc+1])
        line1.set_ydata(u_int[1:nc+2])
        line2.set_ydata(ue)
        plt.title('nc='+str(nc)+', CFL='+str(cfl)+', time ='+str(np.round(t,3)))
        plt.draw(); plt.pause(0.1)

# save exact solution at final time
xe = np.linspace(xmin, xmax, 10000)
ue = uexact(xe,t, uinit)
fname = 'exact.txt'
np.savetxt(fname, np.column_stack([xe, ue]))
print('Saved file ', fname)

# save final time cell-centre solution to a file
fname = 'sol.txt'
np.savetxt(fname, np.column_stack([x, u[1:nc+1]]))
print('Saved file ', fname)

# save final time interface solution to a file
fname = 'sol_int.txt'
np.savetxt(fname, np.column_stack([x_int, u_int[1:nc+2]]))
print('Saved file ', fname)


if args.compute_error == 'yes':
    er1, er2 = compute_error_avg(u[1:nc+1],x, t)         #for error of cell average solution  
    er3, er4 = compute_error(u_int[1:nc+2],x_int, t)     #for error of interface solution
    print('h, L1 error norm, L2 error norm= ')
    print(h, er1, er2, er3, er4)
if args.plot_freq >0:
    plt.show()



