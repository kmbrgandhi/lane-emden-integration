from full_shooting import fitting_method, runge_ketta_shooting_direct, fy, fz

"""
Test 1: Runge Kutta and Shooting Method on Lane Emden equation
Compare printed values to those in the table on page 340 of Stellar Interiors
"""

# Use our methods to find the surface and surface derivatives for different values of n!
lst_of_n = [0, 1, 1.5, 2, 3, 4]
for i in lst_of_n:
	print("n: ", i)
	runge_ketta_shooting_direct(i, fy, fz)
	fitting_method(i, fy, fz)


