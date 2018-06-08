from full_shooting import rungeKetta, get_information_polytrope, solve_equation
from interpolating import evaluate_background
 

def dcdr(r, zeta, delP):
 	if r==0:
 		return 0.1
 	else:
 		return (-1./r) * (3* zeta + (5/3) * delP)  # INSERT GAMMA_1 HERE IN GENERAL

def createDelPdr(P, rho, r_background, M, G, sigma):
	def delPdr(r, zeta, delP, n):
		if r==0:
			return 0.1
		first_part = (1/evaluate_background(r_background, M, r)) * (1/r**(2)) * evaluate_background(r_background, rho, r) * evaluate_background(r_background, M, r) * G
		second_part = (4 * zeta + sigma**(2) * r**(3)/(G * evaluate_background(r_background, M, r)) + delP)
		return -first_part * second_part
	return delPdr


def fitting_adiabatic_radial():
	return

def one_direction_shot_adiabatic_radial(n=2, rho_c = 10**(5), K = 10**(5)):
 	lei= get_information_polytrope(n, rho_c, K)
 	lst_M = lei[-1]
 	lst_P = lei[-2]
 	lst_rho = lei[-3]
 	lst_R = lei[-4]
 	total_mass = lei[-5]
 	surface = lei[-6]
 	G = lei[-7]
 	sigma_guess = 3.9 * (G * total_mass * surface**(-3))**(0.5)
 	zeta_guess = 0.5
 	num_steps = 300
 	hs = surface/(num_steps)
 	for i in range(1000):
 		print(i)
 		R = 0
 		zeta = zeta_guess
 		delP = -3 * (3./5) * zeta
 		delPdr = createDelPdr(lst_P, lst_rho, lst_R, lst_M, G, sigma_guess)
	 	for j in range(num_steps):
	 		R, zeta, delP = rungeKetta(R, zeta, delP, dcdr, delPdr, hs, n)
	 	initial_Y = 4 * zeta + zeta * sigma_guess**(2) * surface**(3)/(G * total_mass) + delP
	 	initial_Z = zeta - 1

	 	delta = 0.001
	 	shiftzeta_guess = zeta_guess*(1 + delta)

	 	zeta = shiftzeta_guess
	 	delP = -3 * (3./5) * zeta
	 	R = 0
	 	for j in range(num_steps):
	 		R, zeta, delP = rungeKetta(R, zeta_guess, delP, dcdr, delPdr, hs, n)
	 	change_Yz = 4*zeta + zeta * sigma_guess**(2) * surface**(3)/(G * total_mass) + delP
	 	change_Zz = zeta - 1

	 	shiftsigma_guess = sigma_guess * (1 + delta)
	 	shiftdelPdr = createDelPdr(lst_P, lst_rho, lst_R, lst_M, G, shiftsigma_guess)
	 	zeta = zeta_guess
	 	delP = -3 * (3./5) * zeta
	 	R = 0
	 	for j in range(num_steps):
	 		R, zeta, delP = rungeKetta(R, zeta_guess, delP, dcdr, delPdr, hs, n)
	 	change_Ysig = 4*zeta + zeta * sigma_guess**(2) * surface**(3)/(G * total_mass) + delP
	 	change_Zsig = zeta - 1

	 	change_zeta, change_sig = solve_equation(-initial_Y, -initial_Z, (change_Yz - initial_Y)/(delta*zeta_guess), (change_Zz - initial_Z)/(delta*zeta_guess), 
			(change_Ysig - initial_Y)/(delta*sigma_guess), (change_Zsig - initial_Z)/(delta*sigma_guess))
	 	if(abs(change_sig) < 10**(-7) and abs(change_zeta) < 10**(-7)):
	 		break
	 	print(change_zeta)
	 	print(change_sig)
	 	zeta_guess = change_zeta + zeta_guess
	 	sigma_guess = change_sig + sigma_guess

one_direction_shot_adiabatic_radial()







 
