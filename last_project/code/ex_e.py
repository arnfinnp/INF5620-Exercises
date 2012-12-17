from dolfin import *
import sys, numpy, time

degree = int(sys.argv[1])
divisions = [int(arg) for arg in sys.argv[2:]]
d = len(divisions)
domain_type = [UnitInterval, UnitSquare, UnitCube]
mesh = domain_type[d-1](*divisions)

V = FunctionSpace(mesh, 'Lagrange', degree)

# Initial condition
u_0 = Constant(0.0)

# Manufactured solution
u_m = Expression('t*x[0]*x[0]*(1.0/2.0 - x[0]/3.0)', t =0.0)

# Set up time steps

t_stop = 1
dt = (1.0/float(sys.argv[2]))**2
t = dt
print 'dx = ', 1.0/float(sys.argv[2]), ', dt = ', dt

# alpha
def alpha(u):
    return 1.0 + u**2

# Define variational problem
rho = 1.0
u_k = interpolate(u_0, V)
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('-rho*pow(x[0],3)/3.0 + rho*pow(x[0],2)/2.0\
			   + pow(t,3)*(8.0*pow(x[0],7)/9.0 - 28.0*pow(x[0],6)/9.0 + 7.0*pow(x[0],5)/2.0 - 5.0*pow(x[0],4)/4.0)\
			   + t*(2.0*x[0] - 1.0)', rho = rho, t=0.0)
a = (rho*inner(u, v) + dt*inner(alpha(u_k)*nabla_grad(u), nabla_grad(v)))*dx
L = (rho*inner(u_k, v) + dt*inner(f, v))*dx

A = assemble(a)
b = None  # variable used for memory savings in assemble calls
u = Function(V)

while t <= t_stop:
    u_m.t = t
    f.t = t
    b = assemble(L, tensor=b)
    solve(A, u.vector(), b)
    um = interpolate(u_m, V)
    t += dt
    u_k.assign(u)

plot(um, title='u_m')
plot(u, title='u')
interactive()