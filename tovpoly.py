import numpy as np
import math
from full_shooting import runge_ketta


gamma = 1.84
c = 2.998 * 10**(10)

def dM(r, rho, M, n):
	return 4 * math.pi * r**(2) * rho

def drho(r, rho, M):
	if rho ==0: 
		return 0
	return -10**(7.36) * rho**(1 - gamma) * gamma**(-1) * G* (rho + (P/(c**(2))))* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))


"""
Perform direct runge_ketta shooting for a particular polytropic index, as well as y', z' functions
@param n the polytropic index
@param f1 how to evaluate y'
@param f2 how to evaluate z'
@return nothing
"""
def shooting_direct(central_density, f1, f2):
	small_x= 0.00001
	y, z = central_density, 0
	x = small_x
	approximate_radius = 1000000
	hs = approximate_radius/250
	lst_x = [0]
	lst_y = [central_density]
	lst_z = [0]
	for i in range(10000):
		x, y, z = rungeKetta(x, y, z, f1, f2, hs, n)
		lst_x.append(x)
		lst_y.append(y)
		lst_z.append(z)
		if abs(y) < 2*hs and hs>(approximate_radius/4000):
			hs = hs/2
		if y<0.000001:
			break
	print("Surface at: ", x)
	print("Total masss: ", z)


shooting_direct(10**(14.25), drho, dM)

