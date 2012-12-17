'''
Poisson equation with Dirichlet conditions
-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2.
f = -6.
'''

from dolfin import *

#mesh = UnitSquare(6, 4)
mesh = UnitCube(6, 4, 5 )

V = FunctionSpace(mesh, 'Lagrange', 1)

# Define boundary conditions
u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
# x[0]: x direction
# x[1]: y direction

def u0_boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u0, u0_boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = f*v*dx

# Compute solution
prm = parameters['krylov_solver']
prm['absolute_tolerance'] = 1E-10
prm['relative_tolerance'] = 1E-6
prm['maximum_iterations'] = 1000
set_log_level(PROGRESS)


u = Function(V)

solve(a == L, u, bc, solver_parameters={'linear_solver':'cg', 'preconditioner': 'ilu'})
#info(solver.parameters, True)
# Plot the solution and mesh
plot(u)
plot(mesh)
print mesh

# Dump solution to file in VTK format
file = File('poisson.pvd')
file << u

# Hold plot
interactive()