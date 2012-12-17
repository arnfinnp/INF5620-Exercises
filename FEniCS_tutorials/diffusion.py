from dolfin import *
import numpy 

mesh = UnitSquare(5, 5)
V = FunctionSpace(mesh, 'Lagrange', 1)

alpha = 3
beta = 1.2

u0 = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', alpha=alpha, beta=beta, t=0)

#u0.t = 0

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u0, boundary)

u_1 = interpolate(u0, V)
# or
#u_1 = project(u0, V) # This results in approximative values at the nodes.

# We could also have used L0 and a0

dt = 0.3 # Time step
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)

a = u*v*dx + dt*inner(nabla_grad(u),nabla_grad(v))*dx
L = (u_1 + dt*f)*v*dx

A = assemble(a)

u = Function(V)
T = 2
t = dt
b = None

while t <= T:
	#b = assemble(L) # This creates a new vector b for each time step. Not very memory friendly!
	b = assemble(L, tensor=b) # This reuses the memory space for b.
	u0.t=t
	bc.apply(A, b)
	solve(A, u.vector(), b)
	plot(u)
	# Verification
	u_e = interpolate(u0, V)
	maxdiff = numpy.abs(u_e.vector().array() - u.vector().array()).max()
	print 'Max error, t=%.2f: %-10.3f' % (t, maxdiff)

	t += dt
	u_1.assign(u)



