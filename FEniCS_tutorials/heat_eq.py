from dolfin import *
import numpy, sys, time

# Setting up the mesh
degree = int(sys.argv[1])
D = float(sys.argv[2])
W = D/2.0
divisions = [int(arg) for arg in sys.argv[3:]]
d = len(divisions)
if d == 1:
	mesh = Interval(divisions[0], -D, 0)
elif d == 2:
	mesh = Rectangle(-W/2, -D, W/2, 0, divisions[0], divisions[1])
elif d == 3:
	mesh = Box(-W/2, -W/2, -D, W/2, W/2, 0, divisions[0], divisions[1], divisions[2])

V = FunctionSpace(mesh, 'Lagrange', 1)

# Setting up Dirichlet boundary conditions
T_R = 0
T_A = 1.0
omega = 2*pi

T_0 = Expression('T_R + T_A*sin(omega*t)', T_R=T_R, T_A=T_A, omega=omega,t=0.0)

def surface(x, on_boundary):
	return on_boundary and abs(x[d-1]) < 1E-14

bc = DirichletBC(V, T_0, surface)

'''
class Kappa(Function):
	def eval(self, value, x):
		"""x: spatial point, value[0]: function value."""
		d = len(x) # no of space dimensions.
		material = 0 # 0: outside, 1: inside.
		kappa_0 = 1
		kappa_1 = 1
		if d == 1:
			if -D/2. < x[d-1] < -D/2. + D/4.:
				material = 1
		elif d == 2:
			if -D/2. < x[d-1] < -D/2. + D/4. and - W/4. < x[0] < W/4.:
				material = 1
		elif d == 3:
			if -D/2. < x[d-1] < -D/2. + D/4. and - W/4. < x[0] < W/4. and -W/4. < x[1] < W/4.:
				material = 1
		value[0] = kappa_0 if material == 0 else kappa_1		
'''

kappa_str = {}
kappa_str[1] = 'x[0} > - D/2 && x[0] < -D/2 + D/4 ? kappa_1 :kappa_0'
kappa_str[2] = 'x[0] > - W/4 && x[0] < W/2 && x[1] > -D/2 && x[1] < -D/2 + D/4 ? kappa_1 : kappa_0'
kappa_str[3] = 'x[0] > - W/4 && x[0] < W/2 && x[1] > -W/4 && x[1] < W/2 && x[2] > -D/2 && x[2] < -D/2 + D/4 ? kappa_1 : kappa_0'

kappa = Expression(kappa_str[d], D=D, W=W, kappa_0=kappa_0, kappa_1=kappa_1)

T_prev = interpolate(Constant(T_R), V)

rho = 1
c = 1
period = 2*pi/omega
t_stop = 5*period
dt = period/20
theta = 1

T = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
a = rho*c*T*v*dx + theta*dt*kappa*inner(nabla_grad(T), nabla_grad(v))*dx
L = (rho*c*T_prev*v + dt*f*v - (1-theta)*dt*kappa*inner(nabla_grad(T), nabla_grad(v)))*dx

A = assemble(a)
b = None
T = Function(V)
t = dt

dummy = Expression('T_R - T_A/2.0 + 2*T_A*(x[%g]+D)' % (d-1),
                   T_R=T_R, T_A=T_A, D=D)

# Make all plot commands inctive
import scitools.misc
plot = scitools.misc.DoNothing(silent=True)
# Need initial dummy plot
viz = plot(dummy, axes=True,
           title='Temperature', wireframe=False)
viz.elevate(-65)
#time.sleep(1)
viz.update(T_1)

import scitools.BoxField
start_pt = [0]*d; start_pt[-1] = -D  # start pt for line plot
import scitools.easyviz as ev

def line_plot():
    """Make a line plot of T along the vertical direction."""
    if T.ufl_element().degree() != 1:
        T2 = interpolate(T, FunctionSpace(mesh, 'Lagrange', 1))
    else:
        T2 = T
    T_box = scitools.BoxField.dolfin_function2BoxField(
            T2, mesh, divisions, uniform_mesh=True)
    #T_box = scitools.BoxField.update_from_dolfin_array(
    #        T.vector().array(), T_box)
    coor, Tval, fixed, snapped = \
            T_box.gridline(start_pt, direction=d-1)

    # Use just one ev.plot command, not hold('on') and two ev.plot
    # etc for smooth movie on the screen
    if kappa_0 == kappa_1:  # analytical solution?
        ev.plot(coor, Tval, 'r-',
                coor, T_exact(coor), 'b-',
                axis=[-D, 0, T_R-T_A, T_R+T_A],
                xlabel='depth', ylabel='temperature',
                legend=['numerical', 'exact, const kappa=%g' % kappa_0],
                legend_loc='upper left',
                title='t=%.4f' % t)
    else:
        ev.plot(coor, Tval, 'r-',
                axis=[-D, 0, T_R-T_A, T_R+T_A],
                xlabel='depth', ylabel='temperature',
                title='t=%.4f' % t)

    ev.savefig('tmp_%04d.png' % counter)
    time.sleep(0.1)

while t <= t_stop:
	b = assemble(L, tensor=b)
	t_0.t = t
	bc.apply(A, b)
	solve(A, T.vector(), b)

	viz.update(T)
    #viz.write_ps('diff2_%04d' % counter)
    line_plot()

	t += dt
	T_prev.assign(T)

viz = plot(T, title='Final solution')
viz.elevate(-65)
viz.azimuth(-40)
viz.update(T)

interactive()