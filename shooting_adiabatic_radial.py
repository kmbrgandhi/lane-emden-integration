from full_shooting import rungeKetta, get_information_polytrope, solve_equation
from interpolating import evaluate_background
 

def dcdr(r, zeta, delP):
 	if r==0:
 		return -0.01
 	else:
 		return -(1./r) * (3* zeta + 0.6 * delP)  # INSERT GAMMA_1 HERE IN GENERAL
Pdict = {}
Mdict = {}
rhodict = {}

def createDelPdr(P, rho, r_background, M, G, sigma, surface_real):
	def delPdr(r, zeta, delP, n):
		if r in Pdict:
			P_value = Pdict[r]
		else:
			P_value = evaluate_background(r_background, P, r)
			Pdict[r] = P_value

		if r in Mdict:
			M_value = Mdict[r]
		else:
			M_value = evaluate_background(r_background, M, r)
			Mdict[r] = M_value
		if r in rhodict:
			rho_value = rhodict[r]
		else:
			rho_value = evaluate_background(r_background, rho, r)
			rhodict[r] = rho_value
		if r==0 or (abs(r - surface_real) < 0.00001):
			return -0.01
		first_part = (1/P_value) * (1/r**(2)) * rho_value * M_value * G
		second_part = (4 * zeta + (zeta * sigma**(2) * r**(3)/(G * M_value)) + delP)
		return (first_part * second_part)
	return delPdr

def fitting_method_subpart(zeta, sigma, f1, f2, lst_x, lst_y, lst_z, surface, G, total_mass):
	fitting_point = surface/2.0
	xout = 0
	yout = zeta
	zout = -5 * zeta
	num_steps = 100
	hsout = fitting_point/(num_steps)
	lst_x.append(xout)
	lst_y.append(yout)
	lst_z.append(zout)
	for j in range(num_steps):
		if xout + hsout - fitting_point > -0.00000001:
			hsout = fitting_point - xout
			xout, yout, zout = rungeKetta(xout, yout, zout, f1, f2, hsout, 2)
			lst_x.append(xout)
			lst_y.append(yout)
			lst_z.append(zout)
			break
		else:
			xout, yout, zout = rungeKetta(xout, yout, zout, f1, f2, hsout, 2)
			lst_x.append(xout)
			lst_y.append(yout)
			lst_z.append(zout)
		
	to_reverse_x = []
	to_reverse_y = []
	to_reverse_z = []
	xin = surface
	yin = 1
	zin = -(4 + (sigma**(2) * surface**(3)/(G * total_mass)))
	hsin = -fitting_point/(num_steps)
	to_reverse_x.append(xin)
	to_reverse_y.append(yin)
	to_reverse_z.append(zin)
	for i in range(num_steps):
		if fitting_point - (xin + hsin) > -0.000000001:
			hsin = fitting_point - xin
			xin, yin, zin = rungeKetta(xin, yin, zin, f1, f2, hsin, 2)
			break
		else:
			xin, yin, zin = rungeKetta(xin, yin, zin, f1, f2, hsin, 2)
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

def fitting_adiabatic_radial(n=2, rho_c = 10**(5), K = 10**(5)):
	P_dict = {}
	rho_dict = {}
	M_dict = {}
	lei= get_information_polytrope(n)
	lst_M = lei[-1]
	lst_P = lei[-2]
	lst_rho = lei[-3]
	lst_R = lei[-4]
	total_mass = lei[-5]
	surface = lei[-6]
	G = lei[-7]
	sigma_guess = (3.9 * G * total_mass * surface**(-3))**(0.5)
	zeta_guess = 0.5
	num_steps = 1000
	hs = surface/(num_steps)

	final_R = []
	this_R = []
	final_zeta = []
	this_zeta= []
	final_delP= []
	this_delP = []
	hs = 0.001
	for i in range(300):
		this_R = []
		this_zeta = []
		this_delP = []
		delPdr = createDelPdr(lst_P, lst_rho, lst_R, lst_M, G, sigma_guess, surface)
		initialY, initialZ = fitting_method_subpart(zeta_guess, sigma_guess, dcdr, delPdr, this_R, this_zeta, this_delP, surface, G, total_mass)
		delta= 0.01
		csg = zeta_guess * (1 + delta)
		csd = sigma_guess * (1 + delta)
		cYz, cZz = fitting_method_subpart(csg, sigma_guess, dcdr, delPdr, [], [], [], surface, G, total_mass)
		newdelPdr = createDelPdr(lst_P, lst_rho, lst_R, lst_M, G, csd, surface)
		cYp, cZp = fitting_method_subpart(zeta_guess, csd, dcdr, newdelPdr, [], [], [], surface, G, total_mass)
		change_z, change_s = solve_equation(-initialY, -initialZ, (cYz - initialY)/(delta*zeta_guess), (cZz - initialZ)/(delta*zeta_guess), (cYp - initialY)/(delta*sigma_guess), (cZp - initialZ)/(delta*sigma_guess))
		if(abs(change_z) < 10**(-7) and abs(change_s) < 10**(-7)):
			break
		zeta_guess = zeta_guess + change_z
		sigma_guess = sigma_guess + change_s
		print(zeta_guess)
		print(sigma_guess)
	final_R = this_R
	final_zeta = this_zeta
	final_delP = this_delP

def one_direction_shot_adiabatic_radial(n=2):
 	lei= get_information_polytrope(n=n)
 	lst_M = lei[-1]
 	lst_P = lei[-2]
 	lst_rho = lei[-3]
 	lst_R = lei[-4]
 	total_mass = lei[-5]
 	surface = lei[-6]
 	G = lei[-7]
 	sigma_guess = (3.9 * G * total_mass * surface**(-3))
 	zeta_outside = 1
 	num_steps =300
 	hs = -surface/(num_steps)
 	for i in range(1000):
 		print(i)
 		R = surface
 		zeta = zeta_outside
 		delP = -(((sigma_guess * surface**(3))/(G * total_mass)) * zeta + 4 * zeta)
 		print(delP)
 		delPdr = createDelPdr(lst_P, lst_rho, lst_R, lst_M, G, sigma_guess, surface)
	 	for j in range(num_steps):
	 		R, zeta, delP = rungeKetta(R, zeta, delP, dcdr, delPdr, hs, n)
	 	bc = 3 * zeta + 0.6 * delP
	 	print(R)
	 	print(zeta)
	 	print(delP)
	 	delta = 0.001

	 	shiftsigma_guess = sigma_guess * (1 + delta)
	 	R = surface
	 	zeta = zeta_outside
 		delP = -(((shiftsigma_guess * surface**(3))/(G * total_mass)) * zeta + 4 * zeta)
 		newdelPdr = createDelPdr(lst_P, lst_rho, lst_R, lst_M, G, shiftsigma_guess, surface)
	 	for j in range(num_steps):
	 		R, zeta, delP = rungeKetta(R, zeta, delP, dcdr, newdelPdr, hs, n)
	 	sbc = 3 * zeta + 0.6 * delP

	 	print(bc)
	 	print(sbc)
	 	break

fitting_adiabatic_radial()
#one_direction_shot_adiabatic_radial()







 
