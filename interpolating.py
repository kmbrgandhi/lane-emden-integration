
"""
Evaluate a background mesh at a particular value using linear interpolation
@param indep: an array of floats corresponding to the values of the independent variable at which we know the value of the function: sorted increasing
@param dep: an array of floats corresponding to the values of the function at the corresponding values of independent variable in indep
@param point: the value at which we wish to evaluate our function: 
@return the value of the function at point
"""
def evaluate_background(indep, dep, point):
	i = 0
	found = False;
	for i in range(len(indep)):
		if(abs(indep[i] -point)<0.000001):
			return dep[i]
		elif indep[i] > point:
			found = True;
			break
	if not found:
		return dep[len(indep) - 1]
	elif i ==0:
		return dep[0]
	else:
		part1 = point - indep[i-1]
		part2 = indep[i] - point
		comput1 = part2 * dep[i-1]
		comput2 = part1 * dep[i]
		full_sum = comput1 + comput2
		return full_sum/(part1 + part2)


def test_interpolating():
	background = [0, 0.1, 0.3, 0.6, 0.85, 1]
	values = [2, 4, 3, 7, 2, 1]
	print(evaluate_background(background, values, 0))
	print("should be: ", 2)
	print(evaluate_background(background, values, 0.6))
	print("should be: ", 7)
	print(evaluate_background(background, values, 1))
	print("should be: ", 1)
	print(evaluate_background(background, values, 0.2))
	print("should be: ", 0.5 * 4 + 0.5 * 3)
	print(evaluate_background(background, values, 0.35))
	print("should be: ", (5./6) * 3 + (1./6) * 7)
	print(evaluate_background(background, values, 0.8))
	print("should be: ", 0.8 * 2 + 0.2 * 7)
