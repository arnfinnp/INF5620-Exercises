from dolfin import *
import sys, numpy, time

degree = int(sys.argv[1])
divisions = [int(arg) for arg in sys.argv[2:]]
d = len(divisions)
domain_type = [UnitInterval, UnitSquare, UnitCube]
mesh = domain_type[d-1](*divisions)

V = FunctionSpace(mesh, 'Lagrange', degree)

# Initial condition
sigma = 0.2

u_0 = Expression('exp(-(1.0/(2*sigma*sigma))*(x[0]*x[0] + x[1]*x[1]))', sigma=sigma)

# Set up time steps

dt = (1.0/float(sys.argv[2]))**2
t = dt
t_stop = 0.8
print 'dx = ', 1.0/float(sys.argv[2]), ', dt = ', dt

# alpha
def alpha(u):
    beta = 5.0
    return 1.0 + beta*u**2

# Define variational problem
rho = 1.0
u_k = interpolate(u_0, V)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)
a = (rho*inner(u, v) + dt*inner(alpha(u_k)*nabla_grad(u), nabla_grad(v)))*dx
L = (rho*inner(u_k, v) + dt*inner(f, v))*dx

A = assemble(a)
b = None  # variable used for memory savings in assemble calls
u = Function(V)

while t <= t_stop:
    b = assemble(L, tensor=b)
    solve(A, u.vector(), b)
    plot(u, wireframe = True, scale = 1.0)
    t += dt
    u_k.assign(u)
    print t

#plot(u, title='u')
interactive()