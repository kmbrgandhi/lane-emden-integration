
"""
Evaluate a background mesh at a particular value using linear interpolation
@param indep: an array of floats corresponding to the values of the independent variable at which we know the value of the function: sorted increasing
@param dep: an array of floats corresponding to the values of the function at the corresponding values of independent variable in indep
@param point: the value at which we wish to evaluate our function: 
@return the value of the function at point
"""
def evaluate_background(indep, dep, point):
	i = 0
	for i in range(len(indep)):
		if(abs(indep[i] -point)<0.000001):
			return dep[i]
		elif indep[i] > point:
			break
	if i ==0:
		return dep[0]
	else:
		part1 = point - indep[i-1]
		part2 = indep[i] - point
		comput1 = part1 * dep[i-1]
		comput2 = part2 * dep[i]
		full_sum = comput1 + comput2
		return full_sum/(part1 + part2)