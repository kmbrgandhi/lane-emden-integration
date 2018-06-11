import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as datetime


# Number of grid points
N = 100
dx = 1/(N+1)
x = np.linspace(dx, 1-dx, N)

# Exact solution
f_e = -8*x+10

# Tolerance
tol = 1e-6
diff = 1
it_error = []

# Difference first round
diff0 = 0

# Boundary conditions
f0 = 10
f1 = 2

# Diagonal elements
d = 2

# Off-diagonal elements
e = -1

# Intitialize matrix and set elements
A = np.zeros([N, N])
for i in range(N):
    A[i, i] = d
    if i < N-1:
        A[i, i+1] = e
    if i > 0:
        A[i, i-1] = e
        

# Initialize b-vector
b = np.zeros(N)
b[0] = f0
b[N-1] = f1
        
# Inititializing solution array, use zeros as initial guess
f = np.zeros(N)

n_iter = 0

while diff > tol*diff0:
    f_new = np.zeros(N)
    for i in range(N):
        if i == 0:
            f_new[i] = (b[i] - A[i, i+1]*f[i+1])/A[i, i]
        elif i == N-1:
            f_new[i] = (b[i] - A[i, i-1]*f_new[i-1])/A[i, i]
        else:
            f_new[i] = (b[i] - A[i, i+1]*f[i+1] - A[i, i-1]*f_new[i-1])/A[i, i]


    # Determine 2-norm of differnce between iterations
    diff = np.sqrt(sum((f_new-f)**2))
    
    # Set 2-norm of difference in first iteration
    if diff0 == 0:
        diff0 = diff
        
    it_error.append(np.sqrt(sum((f_new-f_e)**2)))
        
    # Update f
    f = f_new

    n_iter += 1