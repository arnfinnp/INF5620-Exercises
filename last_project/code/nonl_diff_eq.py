from dolfin import *
import sys, numpy, time
import scitools.BoxField

degree = int(sys.argv[1])
divisions = [int(arg) for arg in sys.argv[2:]]
d = len(divisions)
domain_type = [UnitInterval, UnitSquare, UnitCube]
mesh = domain_type[d-1](*divisions)
V = FunctionSpace(mesh, 'Lagrange', degree)

# Initial condition
u_0 = Expression('cos(pi*x[0])')

# Exact solution
u_e = Expression('exp(-pi*pi*t)*cos(pi*x[0])', t=0.0)

# Set up time steps

t_stop = 0.07
dt = (1.0/float(sys.argv[2]))**2
t = dt
print 'dx = ', 1.0/float(sys.argv[2]), ', dt = ', dt

# alpha
def alpha(u):
    return Expression('1')

# Define variational problem
rho = 1
u_k = interpolate(u_0, V)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
a = (rho*inner(u, v) + dt*inner(alpha(u_k)*nabla_grad(u), nabla_grad(v)))*dx
L = (rho*inner(u_k, v) + dt*inner(f, v))*dx

A = assemble(a)
b = None  # variable used for memory savings in assemble calls

u = Function(V)

while t <= t_stop:
    u_e.t = t
    b = assemble(L, tensor=b)
    solve(A, u.vector(), b)
    ue = interpolate(u_e, V)
    if t > (0.05 - dt) and t < (0.05 + dt):
        e = ue.vector().array() - u.vector().array()
        E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)
        print 'error = ', E, 't = ', t, 'E/h = ', E/dt
    t += dt
    u_k.assign(u)