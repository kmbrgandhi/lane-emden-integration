import numpy as np
import math
from full_shooting import rungeKetta
import scipy.optimize
import matplotlib.pyplot as plt
import skyrme_group1
import eos_bps


# all values in CGS here
c = 2.99792458* 10**(10)
G = 6.67408 * 10**(-8)
K = 10**(7.36)
h = 6.62606885 * 10**(-27)
hbar = 1.0545716 * 10**(-27)
K_nr = ((1.004 * 10**(13))/(1.2**(5./3)))
K_r = ((1.243 * 10**(15))/(1.2**(4./3)))
mass_sun = 1.98847 * 10**(33)

# all values are in MeV or associated units (MeV fm)
hbarc = 197.3269631 # MeV
amu_mass = 931.5
mass_neutron = 938.27231
mass_proton = 939.56563
mass_muon_mev = 105.6583715
hbaroverm = hbarc**(2)/(2*amu_mass)


conversion_pressure = 1.6022 * 10**(33)
conversion_density = 1.7827 * 10**(12)
pi_multiplier= (1.5 * math.pi * math.pi)**(2./3)
gamma = 1.84


names, lst_params, lst_central_14, lst_central_max, correct_radii, correct_masses = skyrme_group1.get_skyrme_info()
eos_rho, eos_P = eos_bps.get_eos_manual()

def asymmetry(n):
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	term1 = 1./3 * hbaroverm * pi_multiplier * n**(2./3)
	term2 = -0.125 * t0 * (2 * x0 + 1) * n
	term3 = -(1./24) * pi_multiplier * (3 * t1 * x1 - t2 * (5*x2 + 4))* n**(5./3)
	term4 = -(1./48) * t3*(2*x3 + 1) * n**(alpha + 1)
	return term1 + term2 + term3 + term4

def get_mean_mass(n):
	I = get_I(n)
	x_p = 0.5 * (1-I)
	x_n = 1 - x_p
	sp = skyrme_params
	ms = sp[9]
	mv = sp[10]
	lhs = ((1+I)/(ms**(2))) - (I/(mv**(2)))
	lhs2 = ((1-I)/(ms**(2))) + (I/(mv**(2)))
	effective_neutron_ratio = (1/(lhs))**(0.5)
	effective_proton_ratio = (1/(lhs2))**(0.5)

def get_xp_eq(n):
	def xp_eq(xp):
		value = xp**(1./3) * (hbarc* (3 * math.pi*math.pi * n)**(1./3))
		value2 = 4 * asymmetry(n) * (1 - 2*xp)
		return (value - value2)
	return xp_eq

def get_xmu_eq(n):
	def xmu_eq(xmu):
		internalterm = (mass_muon_mev/hbarc)**(2) * 1./((3 * math.pi * math.pi *n)**(2./3))
		internalterm = internalterm + xmu**(2./3)
		lhs = (1 - 2*xmu - 2 * (internalterm)**(1.5)) * 4 * asymmetry(n)
		rhsinternal = mass_muon_mev**(2) + hbarc**(2) * (3 * math.pi * math.pi *n)**(2./3) * xmu**(2./3)
		rhs = rhsinternal**(0.5)
		return rhs - lhs
	return xmu_eq

def use_neutrality_simple(n):
	xp_function = get_xp_eq(n)
	value = scipy.optimize.brentq(xp_function, 0, 1)
	return value


def use_neutrality_complex(n):
	xmu_function = get_xmu_eq(n)
	value = scipy.optimize.brentq(xmu_function, 0, 1)
	return value

def get_I_simple(n):
	x_p = use_neutrality_simple(n)
	return 1 - 2*x_p

def get_I(n, sc = False, no_errors = False):

	try:
		x_mu = use_neutrality_complex(n)
		term1 = (mass_muon_mev/hbarc)**(2) * 1./((3 * math.pi * math.pi * n)**(2./3))
		x_e = (term1 + x_mu**(2./3))**(1.5)
		if hbarc * (3 * math.pi * math.pi *n)**(1./3) * x_e**(1./3) < mass_muon_mev:
			return get_I_simple(n)
		x_p = x_e + x_mu
		x_n = 1 - x_mu - x_e - x_p
		if sc:
			return x_p, x_e, x_mu, x_n
		return 1 - 2*x_p
	except Exception as e:
		return get_I_simple(n)
def E(n):
	I = get_I(n)
	"""
	x_p = 0.5 * (1-I)
	x_n = 1 - x_p
	real_mean_mass = x_p * mass_proton + x_n * mass_neutron
	"""
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]

	term1 = 0.6 * hbaroverm*pi_multiplier * n**(2./3) * F(5./3, I)
	term2 = 0.125 * t0 * n * (2*(x0 + 2) - (2*x0 + 1)*F(2, I))
	term3 = (1./48) * t3 * n**(alpha + 1) * (2*(x3 + 2) - (2*x3 + 1)*F(2, I))
	multiplicand = (t1*(x1 + 2) + t2*(x2 + 2))*F(5./3, I) + 0.5*(t2*(2*x2 + 1)-t1*(2*x1 + 1)) * F(8./3, I)
	term4 = (3./40) * pi_multiplier * n**(5./3) * multiplicand
	return term1 + term2 + term3 + term4

def P(n, no_errors = False):
	I = get_I(n, no_errors = no_errors)
	"""
	x_p = 0.5 * (1-I)
	x_n = 1 - x_p
	n_p = n*x_p
	n_n = n* x_n
	real_mean_mass = x_p * mass_proton + x_n * mass_neutron
	"""
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	term1 = 0.4 * hbaroverm*pi_multiplier * n**(5./3) * F(5./3, I)
	term2 = 0.125 * t0 * n**(2) * (2*(x0 + 2) - (2*x0 + 1)*F(2, I))
	term3 = ((alpha + 1)/48.) * t3 * n**(alpha + 2) * (2*(x3 + 2) - (2*x3 + 1)*F(2, I))
	multiplicand = (t1*(x1 + 2) + t2*(x2 + 2))*F(5./3, I) + 0.5*(t2*(2*x2 + 1)-t1*(2*x1 + 1)) * F(8./3, I)
	term4 = (1./8) * pi_multiplier * n**(8./3) * multiplicand
	return term1 + term2 + term3 + term4

def dPdn(n):
	I = get_I(n)
	"""
	x_p = 0.5 * (1-I)
	x_n = 1 - x_p
	n_p = n*x_p
	n_n = n* x_n
	real_mean_mass = x_p * mass_proton + x_n * mass_neutron
	"""
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	term1 = 2./3 * hbaroverm*pi_multiplier * n**(2./3) * F(5./3, I)
	term2 = 0.25 * t0 * n * (2*(x0 + 2) - (2*x0 + 1)*F(2, I))
	term3 = (((alpha + 1)*(alpha+2))/48.) * t3 * n**(alpha + 1) * (2*(x3 + 2) - (2*x3 + 1)*F(2, I))
	multiplicand = (t1*(x1 + 2) + t2*(x2 + 2))*F(5./3, I) + 0.5*(t2*(2*x2 + 1)-t1*(2*x1 + 1)) * F(8./3, I)
	term4 = (1./3) * pi_multiplier * n**(5./3) * multiplicand
	return term1 + term2 + term3 + term4

def get_v_adiabatic(n):
	if n>0:
		top_term = dPdn(n)
		I = get_I(n)
		x_p = 0.5 * (1-I)
		x_n = 1 - x_p
		bot_term = x_p * mass_proton + x_n * mass_neutron + E(n) + (P(n)/n)
		ratio = (top_term/bot_term)**(0.5) 
		return ratio * c
	else:
		return 0

def rho(n):
	return conversion_density * energy_density(n)

def energy_density(n):
	I = get_I(n)
	x_p = 0.5 * (1-I)
	x_n = 1 - x_p
	n_p = n*x_p
	n_n = n* x_n
	energy_density = n*E(n) + n_p * mass_proton + n_n * mass_neutron
	return energy_density

def extract_n_pressure(Pval):
	def desired_fn(n):
		return P(n, no_errors = True) - Pval
	value = scipy.optimize.brentq(desired_fn, 0, 10)
	return value

def getrho(P):
	P_to_check = P/(conversion_pressure)
	value_n = extract_n_pressure(P_to_check)
	rho_val = rho(value_n)
	return rho_val

def dM(r, rho, M, n):
	return 4 * math.pi * r**(2) * rho
def getPpoly(rho, gamma):
	return K * abs(rho)**(gamma)


def getRhoBSM(P, interp = 1.5):
	if (P-eos_P[0] < 0.01) or (P - eos_P[-1] > -0.01):
		raise ValueError('We are out of range of the BSM eos.')
	else:
		for i in range(len(eos_P)):
			if abs(eos_P[i] - P) < 0.1:
				return eos_P[i]

		first_greater = 0
		for i in range(len(eos_P)):
			if eos_P[i] > P:
				first_greater = i
				break

		index_1 = first_greater - 1
		index_2 = first_greater
		try:
			# polynomial fitting
			a_coeff = (eos_P[index_2] - eos_P[index_1])/(eos_rho[index_1]**(interp) - eos_rho[index_2]**(interp))
			b_coeff = 0.5 * (eos_P[index_1] - a_coeff * eos_rho[index_1]**(interp) + eos_P[index_2] - a_coeff * eos_rho[index_2]**(interp))
			return ((P - b_coeff)/a_coeff)**(1./interp)
		except:
			# linear interpolation
			diff = eos_P[index_2] - eos_P[index_1]
			weight_1 = (eos_P[index_2] - P)/diff
			weight_2 = (P - eos_P[index_1])/diff
			return ((weight_1 * eos_rho[index_1]) + (weight_2 * eos_rho[index_2]))





def getRhoPoly(P):
	if P<0:
		return 0
	else:
		K_actual = 5.38 * 10**(9)
		gamma_actual = 5./3
		return (K_actual**(-1) * P)**(1./gamma_actual)

def getRhoPolyextra(P):
	if P<0:
		return 0
	else:
		#include more cases here to cover crust if need be
		value_rel = (K_r**(-1) * P)**(3./4)
		if value_rel < 10**(5):
			return (K_nr**(-1) * P)**(3./5)
		else:
			return value_rel
		#return (K**(-1) * P)**(1./gamma)

def drhoPoly(r, rho, M):
	if rho <=0: 
		return 0
	value = -K**(-1) * abs(rho)**(1 - gamma) * gamma**(-1) * G* (rho + (getPpoly(rho, gamma)/(c**(2))))* ((M + 4 * math.pi * r**(3) * (getPpoly(rho, gamma)/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))
	return value

def dMSkyrme(r, P, M, n):
	surface_area = 4 * math.pi * r**(2)

	try:
		value = getrho(P)
		if getrho(P) > 3*10**(14):
			return surface_area * getrho(P)
		else:
			try:
				return surface_area* getRhoBSM(P)
			except: 
				return surface_area * getRhoPoly(P)
	except:
		try:
			return surface_area * getRhoBSM(P)
		except:
			return surface_area * getRhoPoly(P)

def dPSkyrme(r, P, M):
	if P==0:
		return -2.0*10**(29)
	try:
		if getrho(P) > 3*10**(14):
			return -G* (P/c**(2) + getrho(P))* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))
		elif getrho(P)>10**(14):
			weight = (getrho(P) - 10**(14))/(2 * 10**(14))
			emphasize_bsm_weight = 0.8 * weight
			mean_density = emphasize_bsm_weight * getrho(P) + (1-emphasize_bsm_weight) * getRhoBSM(P)
			return -G* (P/c**(2) + mean_density)* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))
		else:
			try:
				return -G* (P/c**(2) + getRhoBSM(P))* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))
			except: 
				return -G* (P/c**(2) + getRhoPoly(P))* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))
	except Exception as e:
		try:
			return -G* (P/c**(2) + getRhoBSM(P))* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))
		except:
			return -G* (P/c**(2) + getRhoPoly(P))* ((M + 4 * math.pi * r**(3) * (P/(c**(2))))/(r**(2) * (1 - ((2*G*M)/(c**(2) * r)))))

def F(m, I):
	return (0.5 * (1 + I)**(m) + 0.5 * (1-I)**(m))

def get_n_copies(str, n):
	if n==1:
		return str
	else:
		return str + get_n_copies(str, n-1)

def shooting_skyrme(central_pressure):
	small_x = 0.00001
	y, z = central_pressure, 0.0
	x = small_x
	approximate_radius = 2*10**(6)
	hs = approximate_radius/1000
	hs = 1.15988*10**(2)
	lst_x = [0.0]
	lst_y = [central_pressure]
	lst_z = [0.0]
	lst_logbrunt = []
	lst_brunt = []
	lst_x_logbrunt= []
	lst_x_brunt = [] 
	lst_speed_sound = []
	lst_eq_speed = []
	lst_n = []
	lst_xp = []
	lst_xe = []
	lst_xmu = []
	lst_xn = []
	lst_g = []
	lst_rho2 = []
	for i in range(10000):
		x, y, z = rungeKetta(x, y, z, dPSkyrme, dMSkyrme, hs, gamma)
		print(x)
		lst_x.append(x)
		lst_y.append(y)
		lst_z.append(z)
		g = (z * G)/(x**(2)) 
		lst_g.append(g**(2))
		try:
			curr_n = extract_n_pressure(y/conversion_pressure)
			if rho(curr_n) > 1.6*10**(14):
				xp, xe, xmu, xn = get_I(curr_n, sc = True)
				lst_xp.append(xp*curr_n*10**(39))
				lst_xe.append(xe*curr_n*10**(39))
				lst_xmu.append(xmu*curr_n*10**(39))
				lst_xn.append(xn*curr_n*10**(39))
				lst_n.append(curr_n*10**(39))

				adiabatic_speed = get_v_adiabatic(curr_n)
				xdiff, ydiff, zdiff = rungeKetta(x, y, z, dPSkyrme, dMSkyrme, 10, gamma)
				n_diff = extract_n_pressure(ydiff/conversion_pressure)
				equilibrium_ratio =  (ydiff - y)/(conversion_pressure*(energy_density(n_diff) - energy_density(curr_n)))
				if adiabatic_speed/c > 1:
					print(adiabatic_speed/c)
					print('Adiabatic speed suspicious')
				elif equilibrium_ratio > 1:
					print(equilibrium_ratio)
					print('Equilibrium speed suspicious')
				equilibrium_speed = c * equilibrium_ratio
				lst_eq_speed.append(equilibrium_speed**(2))
				brunt = (g**(2) * ((1/equilibrium_speed)**(2) - (1/adiabatic_speed)**(2)))**(0.5)
				lst_x_logbrunt.append(x)
				lst_logbrunt.append(math.log(brunt, 10))
				lst_speed_sound.append(adiabatic_speed**(2))
				lst_x_brunt.append(x)
				lst_brunt.append(brunt**(2))
				lst_rho2.append(rho(curr_n))


		except Exception as e:
			value = 2
		if abs(y) < 2*hs and hs>(approximate_radius/10,000,000):
			hs = hs/2
		if y<0.000001:
			break
	#print(y)

	print("Surface at: ", x/(10**(5)))
	print("Total mass: ", z/mass_sun)
	#return lst_n, lst_xp, lst_xe, lst_xmu, lst_xn
	#return [x/(10**(5)), z/mass_sun]
	lst_mass = [lst_z[i]/mass_sun for i in range(len(lst_z))]

	surface = x
	total_mass = z
	lst_r = lst_x
	lst_m = lst_z
	lst_moverM = [-29.0]
	for i in range(1, len(lst_m)):
		lst_moverM.append(math.log(lst_m[i]/total_mass))
	lst_P = lst_y
	lst_rho = []
	for i in range(len(lst_P)):
		try:
			val = getrho(lst_P[i])
			if val > 3*10**(14):
				lst_rho.append(val)
			else:
				lst_rho.append(getRhoBSM(lst_P[i]))
		except:
			try:
				val = getRhoBSM(lst_P[i])
				lst_rho.append(val)
			except:
				val = getRhoPoly(lst_P[i])
				lst_rho.append(val)
	adiabatic_exp = []
	adiabatic_val = 0
	for i in range(len(lst_rho)):
		y = lst_P[i]
		try:
			curr_n = extract_n_pressure(y/conversion_pressure)
			adiabatic_speed = get_v_adiabatic(curr_n)
			adiabatic_val = adiabatic_speed * (lst_rho[i]/y)
			adiabatic_exp.append(adiabatic_val)
		except:
			if i!=(len(lst_rho)-1):
				top = lst_P[i+1] - lst_P[i-1]
				bot = lst_rho[i+1] - lst_rho[i-1]
				adiabatic_val = top/bot * (lst_rho[i]/y)
				adiabatic>exp.append(adiabatic_val)
			else:
				adiabatic_exp.append(adiabatic_val)
	deriv_diff = []
	diff_val = 0
	for i in range(len(lst_r)):
		if i!=(len(lst_r)-1):
			top1 = lst_P[i+1] - lst_P[i]
			top2 = lst_rho[i+1] - lst_rho[i]
			bot = lst_r[i+1] - lst_r[i]
			term1 = (1./adiabatic_exp[i]) * (top1/bot) * (lst_r[i]/lst_P[i])
			term2 = (top2/bot) * (lst_r[i]/lst_rho[i])
			diff_val = term1 - term2
			deriv_diff.append(diff_val)
		else:
			deriv_diff.append(diff_val)

	lst_Rminusr = [surface - lst_x[i] for i in range(len(lst_r))]
	print(len(lst_r))
	print(len(lst_moverM))
	print(len(lst_P))
	print(len(lst_rho))
	print(len(adiabatic_exp))
	print(len(deriv_diff))
	print(len(lst_Rminusr))
	print('writing')
	file = open('skyrme.fgong','w')
	nl = "\n"
	file.write("fee" + nl)
	file.write("fi" + nl)
	file.write("fo" + nl)
	file.write("fum" + nl)
	separator = " "
	file.write(str(len(lst_r)) + separator + "15" + separator + "40" + separator + "1300" + nl)
	empty_val = str(0.0) + separator
	file.write(str(total_mass) + separator + str(surface) + separator + get_n_copies(empty_val, 12) + str(0.0) + nl)
	for i in range(len(lst_r)):
		print(i)
		full_string = ""
		full_string = full_string + str(lst_r[i]) + separator + str(lst_moverM[i]) + separator
		full_string = full_string + empty_val
		full_string = full_string + str(lst_P[i]) + separator + str(lst_rho[i]) + separator
		full_string = full_string + get_n_copies(empty_val, 4)
		full_string = full_string + str(adiabatic_exp[i]) + separator
		full_string = full_string + get_n_copies(empty_val, 4)
		full_string = full_string + str(deriv_diff[i]) + separator
		full_string = full_string + get_n_copies(empty_val, 2)
		full_string = full_string + str(lst_Rminusr[i]) + separator
		full_string = full_string + get_n_copies(empty_val, 21)
		full_string = full_string + str(0.0) + nl
		file.write(full_string)
	file.close()
	return surface, total_mass, lst_r, lst_moverM, lst_P, lst_rho, adiabatic_exp, deriv_diff, lst_Rminusr

def save_to_fgong(which_index, central_pressure):
	skyrme_params = lst_params[which_index]
	# need to add more parameters
	x, brunt, x2, logbrunt, x_rho, rho, P= shooting_skyrme(central_pressure)
	# saving part




def solve_equation(Y, Z, dYdx, dZdx, dYdz, dZdz):
	a = np.array([[dYdx, dYdz], [dZdx, dZdz]])
	b = np.array([Y, Z])
	result = np.linalg.solve(a, b)
	return result[0], result[1]

def fitting_method_skyrme(central_pressure):
	surface_guess = 10**(5) * 12
	mass_guess = 1.4 * 1.989 * 10**(33)
	f1 = dPSkyrme
	f2 = dMSkyrme

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
		initialY, initialZ = fitting_method_subpart_skyrme(central_pressure, surface_guess, mass_guess, f1, f2, this_x, this_y, this_z)
		delta= 0.001
		csg = surface_guess * (1 + delta)
		csd = mass_guess * (1 + delta)
		cxy, cxz = fitting_method_subpart_skyrme(central_pressure, csg, mass_guess, f1, f2, [], [], [])
		czy, czz = fitting_method_subpart_skyrme(central_pressure, surface_guess, csd, f1, f2, [], [], [])
		change_x, change_z = solve_equation(-initialY, -initialZ, (cxy - initialY)/(delta*surface_guess), (cxz - initialZ)/(delta*surface_guess), (czy - initialY)/(delta*mass_guess), (czz - initialZ)/(delta*mass_guess))
		if(abs(change_x) < 5*10**(4) and abs(change_z) < 10**(31)):
			break
		surface_guess = surface_guess + change_x
		mass_guess = mass_guess + change_z
		print(change_x)
		print(change_z)
	final_x = this_x
	final_y = this_y
	final_z = this_z
	print("Surface at: ", surface_guess/(10**(5)))
	print("Total mass: ", mass_guess/mass_sun)
	#return final_x, final_y, final_z, surface_guess, mass_guess
	return [surface/10**(5), mass_guess/mass_sun]

def fitting_method_subpart_skyrme(central_pressure, surface_guess, mass_guess, f1, f2, lst_x, lst_y, lst_z):
	fitting_point = surface_guess/2.0
	small_x = 0.001
	yout, zout = central_pressure, 0
	xout = small_x
	num_steps = 400
	hsout = fitting_point/float(num_steps)
	lst_x.append(0)
	lst_y.append(yout)
	lst_z.append(zout)

	for j in range(num_steps):
		if xout + hsout - fitting_point > -0.00000001:
			hsout = fitting_point - xout
			xout, yout, zout = rungeKetta(xout, yout, zout, f1, f2, hsout, 0)
			lst_x.append(xout)
			lst_y.append(yout)
			lst_z.append(zout)
			break
		else:
			xout, yout, zout = rungeKetta(xout, yout, zout, f1, f2, hsout, 0)
			lst_x.append(xout)
			lst_y.append(yout)
			lst_z.append(zout)

	to_reverse_x = []
	to_reverse_y = []
	to_reverse_z = []
	xin = surface_guess
	yin, zin = 0, mass_guess
	hsin = -fitting_point/float(num_steps)
	to_reverse_x.append(xin)
	to_reverse_y.append(yin)
	to_reverse_z.append(zin)
	for i in range(num_steps):
		if fitting_point - (xin + hsin) > -0.000000001:
			hsin = fitting_point - xin
			xin, yin, zin = rungeKetta(xin, yin, zin, f1, f2, hsin, 0)
			break
		else:
			xin, yin, zin = rungeKetta(xin, yin, zin, f1, f2, hsin, 0)
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

skyrme_params = lst_params[8]
shooting_skyrme(1.36*10**(35))
"""
plt.plot(rho, P)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Density')
plt.ylabel('Pressure')
plt.savefig('testing_rho' + names[8] + ".png")
plt.clf()
"""

#print(asymmetry(0.158))
#print(get_v_adiabatic(0.158)/c)
#print(get_mean_mass(0.158))
#shooting_skyrme(P(lst_central_max[0]) * conversion_pressure)
#fitting_method_skyrme(P(lst_central_14[0]) * conversion_pressure)
#lst_n, lst_xp, lst_xe, lst_xmu, lst_xn = shooting_skyrme(P(lst_central_max[0]) * conversion_pressure)
##print("I'm here")
#plt.plot(lst_n, lst_xmu)
#plt.yscale('log')
#plt.show()
#fitting_method_skyrme(P(lst_central_max[0]) * conversion_pressure)

"""
for i in range(len(lst_params)):
	print(names[i])
	skyrme_params = lst_params[i]
	x, brunt, x2, logbrunt = shooting_skyrme(P(lst_central_14[i]) * conversion_pressure)
	plt.plot(x, brunt, color = "dodgerblue")
	plt.xlabel('Radius (cm)')
	plt.ylabel('Brunt (s^-1)')
	plt.savefig("core_" + names[i] + ".png")
	plt.clf()
	plt.plot(x2, logbrunt)
	plt.xlabel('Radius (cm)')
	plt.ylabel('log Brunt (s^-1)')
	plt.savefig('core_log' + names[i] + ".png")
	plt.clf()




radii = []
max_mass = []
for i in range(len(lst_params)):
	if i%5==0:
		print(names[i])
	skyrme_params = lst_params[i]
	sp = skyrme_params
	x0, x1, x2, x3 = sp[0], sp[1], sp[2], sp[3]
	t0, t1, t2, t3, alpha = sp[4], sp[5], sp[6], sp[7], sp[8]
	result1 = shooting_skyrme(P(lst_central_14[i]) * conversion_pressure)
	surface1 = result1[0]
	result2 = shooting_skyrme(P(lst_central_max[i]) * conversion_pressure)
	mass2 = result2[1]
	radii.append(surface1)
	max_mass.append(mass2)

diff_radii = [round(-(radii[i] - correct_radii[i]), 2) for i in range(len(radii))]
diff_mass = [round(max_mass[i] - correct_masses[i], 2) for i in range(len(radii))]
print(diff_radii)
print(diff_mass)
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
	print("Surface at: ", x/(10**(5)))
	print("Total mass: ", z/mass_sun)

lst = [10**(14.25 + n*0.25) for n in range(0, 9)]

#for i in lst:
#	shooting_direct(i, drhoPoly, dM)

