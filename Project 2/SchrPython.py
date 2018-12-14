import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("eigenvaluesSchr.dat", unpack=True)

eigenAnalyt = data[1]
eigenNum = np.sort(data[0])
globalError = np.linalg.norm(eigenAnalyt[:10] - eigenNum[:10])
for i in range(len(eigenAnalyt)):
#    print i, abs(eigenAnalyt[i] - eigenNum)
    print eigenNum[i], eigenAnalyt[i]

print "Global Error for ten lowest energy levels"
print globalError


n_arr = np.zeros(18)
for i in range(1, 16):
    n_arr[i] = 10*i
n_arr[-1] = 250
n_arr[-2] = 200

globerror = np.array([45.8,7.0,4.9,2.5,1.6,1.1,0.8,0.6,0.5,\
0.4,0.3,0.3,0.2,0.2,0.2,0.1,0.1,0.1])
globerror2 = np.array([6.4,2.4,4.1,4.7,4.9,5.1,5.3,5.3,5.4,5.5])
n_arr2 = np.array([10,30,50,70,90,110,130,150,200,300])

plt.figure()
plt.plot(n_arr, globerror, label=r"$\rho_N=10$")
plt.plot(n_arr2, globerror2, label=r"$\rho_N=5$")
plt.xlabel("N")
plt.ylabel(r"Error $\epsilon_E$")
plt.grid(True)
plt.legend()
plt.show()
