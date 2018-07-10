import numpy as np
import math
from full_shooting import rungeKetta


gamma = 1.84
c = 2.998 * 10**(10)
G = 6.674 * 10**(-8)
K = 10**(7.36)
h = 6.62606885 * 10**(-27)
hbar = 1.0545716 * 10**(-27)
skyrme_params = [-0.95, -5.78, -1.29, -1.56, -1913.6, 439.8, 2697.6, 10592, 0.25]
# this is SkI1

mass_neutron = 1.674929 * 10**(-24)
mass_proton = 1.672623 * 10**(-24)
mass_muon = 1.8835327 * 10**(-25)

def use_neutrality_complex(n):
	return 0

def use_neutrality_simple(n):
	return 0


def extract_n_pressure(P):
	return 0

def getrho(P):
	value_n = extract_n_pressure(P)
	return rho(n, skyrme_params)

def dMSkyrme(r, P, M, n):
	return 4 * math.pi * r**(2) * getrho(P)

def dPSkyrme(r, P, M, n):
	return G* (P/c**(2) + getrho(P))* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))

def dM(r, rho, M, n):
	return 4 * math.pi * r**(2) * rho

def getPpoly(rho, gamma):
	return K * abs(rho)**(gamma)

def drhoPoly(r, rho, M):
	if rho <=0: 
		return 0
	value = -K**(-1) * abs(rho)**(1 - gamma) * gamma**(-1) * G* (rho + (getPpoly(rho, gamma)/(c**(2))))* ((M + 4 * math.pi * r**(3) * (getPpoly(rho, gamma)/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))
	return value

def asymmetry(n, skyrme_params):
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	term1 = 1/3 * hbar**(2)/(2 * mass_neutron) * (1.5 * math.pi * math.pi)**(2./3) * n**(2./3)
	term2 = 0.125 * t0 * (2 * x0 + 1) * n
	term3 = -(1./24) * (1.5 * math.pi * math.pi)**(2./3) * (3 * t1 * x1 - t2 * (5*x2 + 4))* n**(5./3)
	term4 = -(1./48) * t3*(2*x3 + 1) * n**(alpha + 1)
	return term1 + term2 + term3 + term4

def F(m, I):
	return (0.5 * (1 + I)**(m) + 0.5 * (1-I)**(m))

def get_I_complex(n):
	x_mu = use_neutrality_complex(n)
	term1 = (mass_muon * c/hbar)**(2) * 1/((3 * math.pi * math.pi * n)**(2./3))
	x_e = (term1 + x_mu**(2./3))**(1.5)
	x_p = x_e + x_mu
	return 1 - 2*x_p

def get_I(n):
	x_p = use_neutrality_simple(n)
	return 1 - 2*x_p

def E(n, skyrme_params):
	I = get_I(n)
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	term1 = 0.6 * hbar**(2)/(2 * mass_neutron)*(1.5 * math.pi * math.pi)**(2./3) * n**(2./3) * F(5./3, I)
	term2 = 0.125 * t0 * n * (2*(x0 + 2) - (2*x0 + 1)*F(2, I))
	term3 = (1./48) * t3 * n**(alpha + 1) * (2*(x3 + 2) - (2*x3 + 1)*F(2, I))
	multiplicand = (t1*(x1 + 2) + t2*(x2 + 2))*F(5./3, I) + 0.5*(t2*(2*x2 + 1)-t1*(2*x1 + 1)) * F(8./3, I)
	term4 = (3./40) * (1.5 * math.pi * math.pi)**(2./3) * n**(5./3) * multiplicand
	return term1 + term2 + term3 + term4

def P(n, skyrme_params):
	I = get_I(n)
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	term1 = 0.4 * hbar**(2)/(2 * mass_neutron)*(1.5 * math.pi * math.pi)**(2./3) * n**(5./3) * F(5./3, I)
	term2 = 0.125 * t0 * n**(2) * (2*(x0 + 2) - (2*x0 + 1)*F(2, I))
	term3 = ((alpha + 1)/48) * t3 * n**(alpha + 2) * (2*(x3 + 2) - (2*x3 + 1)*F(2, I))
	multiplicand = (t1*(x1 + 2) + t2*(x2 + 2))*F(5./3, I) + 0.5*(t2*(2*x2 + 1)-t1*(2*x1 + 1)) * F(8./3, I)
	term4 = (1./8) * (1.5 * math.pi * math.pi)**(2./3) * n**(8./3) * multiplicand
	return term1 + term2 + term3 + term4

def dPdn(n, skyrme_params):
	I = get_I(n)
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	term1 = 2./3 * hbar**(2)/(2 * mass_neutron)*(1.5 * math.pi * math.pi)**(2./3) * n**(2./3) * F(5./3, I)
	term2 = 0.25 * t0 * n * (2*(x0 + 2) - (2*x0 + 1)*F(2, I))
	term3 = (((alpha + 1)*(alpha+2))/48) * t3 * n**(alpha + 1) * (2*(x3 + 2) - (2*x3 + 1)*F(2, I))
	multiplicand = (t1*(x1 + 2) + t2*(x2 + 2))*F(5./3, I) + 0.5*(t2*(2*x2 + 1)-t1*(2*x1 + 1)) * F(8./3, I)
	term4 = (1./3) * (1.5 * math.pi * math.pi)**(2./3) * n**(5./3) * multiplicand
	return term1 + term2 + term3 + term4

def rho(n, skyrme_params):
	I = get_I(n)
	n_p = n* 0.5 * (1 - I)
	n_n = n* 0.5 * (1 + I)
	epsilon = E(n, skyrme_params) + c**(2) * (n_n * mass_neutron + n_p * mass_proton)
	return epsilon/c**(2)

def depdn(n, skyrme_params):
	return mass_neutron * c**(2) + E(n, skyrme_params) + P(n, skyrme_params)/n

value = dPdn(0.161, skyrme_params)/(depdn(0.161, skyrme_params))
print(value)

def shooting_skyrme(central_pressure):
	small_x = 0.00001
	y, z = central_pressure, 0
	x = small_x
	approximate_radius = 2*10**(6)
	hs = approximate_radius/1000
	lst_x = [0]
	lst_y = [central_pressure]
	lst_z = [0]
	for i in range(10000):
		x, y, z = rungeKetta(x, y, z, dPSkyrme, dMSkyrme, hs, gamma)
		lst_x.append(x)
		lst_y.append(y)
		lst_z.append(z) 
		if abs(y) < 3*hs and hs>(approximate_radius/10000000):
			hs = hs/2
		if y<0.000001:
			break
	print("Surface at: ", x)
	print("Total mass: ", z)

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
	approximate_radius = 2*10**(6)
	hs = approximate_radius/1000
	lst_x = [0]
	lst_y = [central_density]
	lst_z = [0]
	for i in range(10000):
		x, y, z = rungeKetta(x, y, z, f1, f2, hs, gamma)
		lst_x.append(x)
		lst_y.append(y)
		lst_z.append(z) 
		if abs(y) < 3*hs and hs>(approximate_radius/10000000):
			hs = hs/2
		if y<0.000001:
			break
	print("Surface at: ", x)
	print("Total mass: ", z)



lst = [10**(14.25 + n*0.25) for n in range(0, 9)]

#for i in lst:
#	shooting_direct(i, drhoPoly, dM)

