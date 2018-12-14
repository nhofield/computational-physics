import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("eigenvaluesQdots.dat", unpack=True)
data2 = np.loadtxt("eigenvectorQdots.dat", unpack=True)
eigenvectorGstate = data2[0]
rho_vals = data2[1]
eigenAnalyt = data[1]
eigenNum = np.sort(data[0])

print "Eigenvalues G-state to 10th state"
print eigenNum[:10]


plt.figure()
plt.plot(rho_vals, np.sort(eigenvectorGstate))
plt.xlabel(r"$\rho$", fontsize=15)
plt.ylabel(r"$u_i$", fontsize=15)
plt.show()
