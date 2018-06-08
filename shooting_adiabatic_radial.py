from full_shooting import rungeKetta
from interpolating import evaluate_background
 

 def dcdr(r, zeta, delP):
 	if r==0:
 		return 0.1
 	else:
 		return (-1./r) * (3* zeta + 5 * delP)  # INSERT GAMMA_1 HERE for 5

 def createDelPdr(P, rho, r_background, M, G, sigma):
 	def delPdr(r, zeta, delP):
 		first_part = (1/evaluate_background(r_background, M, r)) * (1/r^2) * evaluate_background(r_background, rho, r) * evaluate_background(r_background, M, r) * G
 		second_part = (4 * zeta + sigma**(2) * r**(3)/(G * evaluate_background(r_background, M, r)) + delP)
 		return -first_part * second_part
 	return delPdr

 def fitting_adiabatic_radial():
 	
