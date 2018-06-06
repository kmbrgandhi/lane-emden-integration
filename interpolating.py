
"""
Evaluate a background mesh at a particular value using linear interpolation
@param indep: an array of floats corresponding to the values of the independent variable at which we know the value of the function: sorted increasing
@param dep: an array of floats corresponding to the values of the function at the corresponding values of independent variable in indep
@param point: the value at which we wish to evaluate our function
@return the value of the function at point
"""
def evaluate_background(indep, dep, point):
	if no