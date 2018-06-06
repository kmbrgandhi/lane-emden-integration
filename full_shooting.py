import numpy as np

def rK3(x, y, z, f1, f2, hs, n):
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

def fy(x, y, z):
	return z

def fz(x, y, z, n = 1.5):
	if x==0:
		return -0.33333333
	else:
		return -(abs(y)**(n)) - 2 * (z/x)

def init_computation(small_x, n):
	y = 1 - 1.6 * small_x**(2) + n/120. * small_x**(4) - (n * (8*n-5)/15120.) * small_x**(6)
	z = -2*1.6 * small_x**(2) + 4*n/120. * small_x**(3) - 6 * (n * (8*n -5)/15120. * small_x**(5))
	return y, z

def lane_emden_shooting_direct(n, f1, f2):
	small_x= 0.001
	y, z = init_computation(small_x, n)
	x = small_x
	hs = 0.025
	for i in range(10000):
		x, y, z = rK3(x, y, z, f1, f2, hs, n)
		if abs(y) < 2*hs and hs>0.000001:
			hs = hs/2
		if y<0.000000001:
			break
	print("Surface at: ", x)
	print("Surface Derivative: ", z)

lst_of_n = [0, 1, 1.5, 2, 3, 4]


def solve_equation(Y, Z, dYdx, dZdx, dYdz, dZdz):
	a = np.array([[dYdx, dYdz], [dZdx, dZdz]])
	b = np.array([Y, Z])
	result = np.linalg.solve(a, b)
	return result[0], result[1]

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
	print("Surface at: ", surface_guess)
	print("Surface Derivative: ", surface_derivative)
	return final_x, final_y, final_z, surface_guess, surface_derivative

def fitting_method_subpart(n, surface_guess, surface_derivative, f1, f2, lst_x, lst_y, lst_z):
	fitting_point = surface_guess/2.0
	small_x = 0.001
	yout, zout = init_computation(small_x, n)
	xout = small_x
	hsout = 0.001
	lst_x.append(xout)
	lst_y.append(yout)
	lst_z.append(zout)
	for j in range(10000):
		if xout + hsout - fitting_point > -0.00000001:
			hsout = fitting_point - xout
			xout, yout, zout = rK3(xout, yout, zout, f1, f2, hsout, n)
			lst_x.append(xout)
			lst_y.append(yout)
			lst_z.append(zout)
			break
		else:
			xout, yout, zout = rK3(xout, yout, zout, f1, f2, hsout, n)
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
	for i in range(10000):
		if fitting_point - (xin + hsin) > -0.000000001:
			hsin = fitting_point - xin
			xin, yin, zin = rK3(xin, yin, zin, f1, f2, hsin, n)
			break
		else:
			xin, yin, zin = rK3(xin, yin, zin, f1, f2, hsin, n)
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


for i in lst_of_n:
	print(i)
	lane_emden_shooting_direct(i, fy, fz)
	fitting_method(i, fy, fz)

lst_x, lst_y, lst_z, surface, surface_derivative = fitting_method(1.0, fy, fz)
pressure_lst = []

