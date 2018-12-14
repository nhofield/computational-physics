import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("eigenvalues.dat", unpack=True)
N = data[1]
eigenvals = np.sort(data[0])

plt.figure()
plt.plot(N, eigenvals)
plt.xlabel("N")
plt.ylabel("Eigenvalues")
plt.show()
