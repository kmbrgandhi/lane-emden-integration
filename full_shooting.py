import numpy as np
import math
import matplotlib.pyplot as plt
from interpolating import evaluate_background
import h5py



"""
Perform Runge ketta shooting
@param x the independent variable
@param y the first dependent variable being solved for
@param z the second dependent variable being solved for
@param f1 how to evaluate y'
@param f2 how to evaluate z'
@param hs the step size
@param n the polytropic index
"""
def rungeKetta(x, y, z, f1, f2, hs, n):
	k1 = f1(x, y, z) * hs
	l1 = f2(x, y, z, n) * hs
	k2 = f1(x + (hs/2), y + (k1/2), z + (l1/2)) * hs
	l2 = f2(x + (hs/2), y + (k1/2), z + (l1/2), n) * hs
	k3 = f1(x + (hs/2), y + (k2/2), z + (l2/2)) * hs
	l3 = f2(x + (hs/2), y + (k2/2), z + (l2/2), n) * hs
	k4 = f1(x + hs, y + k3, z + l3) * hs
	l4 = f2(x + hs, y+ k3, z + l3, n) * hs
	y = y + k1/6 + k2/3 + k3/3 + k4/6
	z = z + l1/6 + l2/3 + l3/3 + l4/6
	x = x + hs
	return x, y, z




"""
A particular fy, corresponding to the lane emden equation
@param x the independent variable, a proxy for radius
@param y the first dependent variable being solved for, a proxy for density
@param z the second dependent variable being solved for, a proxy for y'
@return y'
"""
def fy(x, y, z):
	return z
"""
A particular fz, corresonding to the lane emden equation.
Same specs as fy
@param n: the polytropic index
@return z' at these particular x, y, z and polytropic index
"""
def fz(x, y, z, n = 1.5):
	if x==0:
		return -0.33333333
	else:
		return -(abs(y)**(n)) - 2 * (z/x)


"""
Compute initial values of y, y' at a particular value, using the taylor expansion of the lane emden solution around 0
@param small_x the value of x
@param n the polytropic index
@return y, z for this x, n
"""
def init_computation(small_x, n):
	y = 1 - 1.6 * small_x**(2) + n/120. * small_x**(4) - (n * (8*n-5)/15120.) * small_x**(6)
	z = -2*1.6 * small_x**(2) + 4*n/120. * small_x**(3) - 6 * (n * (8*n -5)/15120. * small_x**(5))
	return y, z




"""
Perform direct runge_ketta shooting for a particular polytropic index, as well as y', z' functions
@param n the polytropic index
@param f1 how to evaluate y'
@param f2 how to evaluate z'
@return nothing
"""
def runge_ketta_shooting_direct(n, f1, f2):
	small_x= 0.001
	y, z = init_computation(small_x, n)
	x = small_x
	hs = 0.025
	for i in range(10000):
		x, y, z = rungeKetta(x, y, z, f1, f2, hs, n)
		if abs(y) < 2*hs and hs>0.000001:
			hs = hs/2
		if y<0.000000001:
			break
	print("Surface at: ", x)
	print("Surface Derivative: ", z)



"""
Solve a linear system of equations
@param Y, Z the RHS of the system
@param dYdx, dYdz the coefficients of ∆x and ∆z in the first equation
@param dZdx, dZz the coefficients of ∆x and ∆z in the second equation
@return ∆x and ∆z that solve the equation
"""
def solve_equation(Y, Z, dYdx, dZdx, dYdz, dZdz):
	a = np.array([[dYdx, dYdz], [dZdx, dZdz]])
	b = np.array([Y, Z])
	result = np.linalg.solve(a, b)
	return result[0], result[1]

"""
Perform two-way fitting for the lane emden equation on a particular polytropic index, as well as y', z' functions
@param n the polytropic index
@param f1 how to evaluate y'
@param f2 how to evaluate z'
@return the final values of x, y, z, as well as the surface and surface derivative, in the final model found
"""
def fitting_method(n, f1, f2):
	surface_guess = 3
	if n==0:
		surface_derivative = -1
	else:
		surface_derivative = -1/(n**(2.5))

	final_x = []
	this_x = []
	final_y = []
	this_y= []
	final_z = []
	this_z = []

	for i in range(300):
		this_x = []
		this_y = []
		this_z = []
		initialY, initialZ = fitting_method_subpart(n, surface_guess, surface_derivative, f1, f2, this_x, this_y, this_z)
		delta= 0.001
		csg = surface_guess * (1 + delta)
		csd = surface_derivative * (1 + delta)
		cxy, cxz = fitting_method_subpart(n, csg, surface_derivative, f1, f2, [], [], [])
		czy, czz = fitting_method_subpart(n, surface_guess, csd, f1, f2, [], [], [])
		change_x, change_z = solve_equation(-initialY, -initialZ, (cxy - initialY)/(delta*surface_guess), (cxz - initialZ)/(delta*surface_guess), (czy - initialY)/(delta*surface_derivative), (czz - initialZ)/(delta*surface_derivative))
		if(abs(change_x) < 10**(-7) and abs(change_z) < 10**(-7)):
			break
		surface_guess = surface_guess + change_x
		surface_derivative = surface_derivative + change_z
	final_x = this_x
	final_y = this_y
	final_z = this_z
	#print("Surface at: ", surface_guess)
	#print("Surface Derivative: ", surface_derivative)
	return final_x, final_y, final_z, surface_guess, surface_derivative

"""
Helper function for the fitting method
@param n the polytropic index
@param surface_guess the current surface guess
@param surface_derivative the current guess for the surface derivative
@param f1, f2: how we evaluate y', z'
@param lst_x, lst_y, lst_z lists to fill with the values of x, y, z that we get on this iteration of two-way shooting
"""
def fitting_method_subpart(n, surface_guess, surface_derivative, f1, f2, lst_x, lst_y, lst_z):
	fitting_point = surface_guess/2.0
	small_x = 0.001
	yout, zout = init_computation(small_x, n)
	xout = small_x
	hsout = 0.001
	lst_x.append(0)
	lst_y.append(yout)
	lst_z.append(zout)
	num_steps = 10000
	for j in range(num_steps):
		if xout + hsout - fitting_point > -0.00000001:
			hsout = fitting_point - xout
			xout, yout, zout = rungeKetta(xout, yout, zout, f1, f2, hsout, n)
			lst_x.append(xout)
			lst_y.append(yout)
			lst_z.append(zout)
			break
		else:
			xout, yout, zout = rungeKetta(xout, yout, zout, f1, f2, hsout, n)
			lst_x.append(xout)
			lst_y.append(yout)
			lst_z.append(zout)
		
	to_reverse_x = []
	to_reverse_y = []
	to_reverse_z = []
	xin = surface_guess
	yin, zin = 0, surface_derivative
	hsin = -0.001
	to_reverse_x.append(xin)
	to_reverse_y.append(yin)
	to_reverse_z.append(zin)
	for i in range(num_steps):
		if fitting_point - (xin + hsin) > -0.000000001:
			hsin = fitting_point - xin
			xin, yin, zin = rungeKetta(xin, yin, zin, f1, f2, hsin, n)
			break
		else:
			xin, yin, zin = rungeKetta(xin, yin, zin, f1, f2, hsin, n)
			to_reverse_x.append(xin)
			to_reverse_y.append(yin)
			to_reverse_z.append(zin)
	to_reverse_x.reverse()
	to_reverse_y.reverse()
	to_reverse_z.reverse()
	lst_x.extend(to_reverse_x)
	lst_y.extend(to_reverse_y)
	lst_z.extend(to_reverse_z)
	Y = yin - yout
	Z = zin - zout
	return Y, Z


def generate_polytrope_file(n, filename, Gamma_1):
	lst_xi, lst_theta, lst_theta_deriv, surface, surface_derivative = fitting_method(n, fy, fz)
	f_lst_xi = np.array(lst_xi)
	f_lst_theta = np.array(lst_theta)
	f_lst_theta_deriv = np.array(lst_theta_deriv)
	f = h5py.File(filename)
	f.create_dataset("Theta", data = f_lst_theta)
	f.create_dataset("dTheta", data = f_lst_theta_deriv)
	f.create_dataset("xi", data = f_lst_xi)
	polytropic_index = np.array([n])
	f.attrs.__setitem__("n", len(lst_theta))
	f.attrs.__setitem__("n_d", 0)
	f.attrs.__setitem__("n_poly", polytropic_index)
	f.attrs.__setitem__("Gamma_1", Gamma_1)

	f.close()



def get_information_polytrope(n = 3.0, totalmass = 1.989 * 10**(33), radius = 6.95508 * 10**(10), G = 6.674 * 10**(-8)):
	lst_xi, lst_theta, lst_theta_deriv, surface, surface_derivative = fitting_method(n, fy, fz)
	r_n = radius/(surface)
	P_c = (1/(4 * math.pi * (n + 1) * surface_derivative**(2)))* G*totalmass**(2) * radius**(-4)
	rho_c = -totalmass * radius**(-3) * (surface/(4 * math.pi * surface_derivative))
	lst_R = [lst_xi[i] * r_n for i in range(len(lst_xi))]
	lst_rho = [lst_theta[i]**(n) * rho_c for i in range(len(lst_theta))]
	Kfirst = ((4 * math.pi)/(surface**(n + 1) * (-surface_derivative)**(n - 1)))**(1./n)
	K = Kfirst * G * (1./(n + 1)) * totalmass**(1 - (1./n)) * radius**(-1 + (3./n))
	lst_P = [P_c * lst_theta[i]**(1+n) for i in range(len(lst_theta))]
	lst_M = [-4 * math.pi * r_n**(3) * rho_c * lst_xi[i]**(2) * lst_theta_deriv[i] for i in range(len(lst_xi))]

	return [P_c, r_n, rho_c, K, G, radius, totalmass, lst_R, lst_rho, lst_P, lst_M]

get_information_polytrope()
"""
info = get_information_polytrope(n=1.0)
lst_R = info[7]
lst_rho = info[8]
lst_P = info[9]
lst_M = info[10]
plt.plot(lst_R, lst_rho)
plt.show()
"""

