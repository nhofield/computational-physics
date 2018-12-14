import numpy as np
import time

N = 8
h = 1./N
hh = h*h
A = np.zeros((N,N), float)
for i in range(N):
    A[i][i] = 2.0/hh
for i in range(N-1):
    A[i][i+1] = A[i+1][i] = -1.0/hh

print "Matrix A"
print A
print "Eigenvalues of A"
now = int(round(time.time() * 1000))
eigenv = np.linalg.eig(A)
print eigenv[0]
whataboutnow = int(round(time.time() * 1000))
print "Runtime =%f ms" % (whataboutnow-now)

with open("eigenvaluespython.dat", "w") as outfile:

    for i in range(N):
        outfile.write("%8.6f \n" % eigenv[0][i])
