import numpy as np
import matplotlib.pyplot as plt
import os

L = 20          #Lattice length
T_init = 2.4      #Initial temperature
T_final = 2.4     #Final temperature
h_T = 1         #Temperature step-size
mc_init = 0     #Initial #MC-cycles
mc_final = 5    #Final #MC-cycles
h_mc = 100      #MC-cycle step-size
up = 'd'        #All spins up or down

programName = "mcmcIsing.cpp"
filename = "resultsL"

F = 15          #Fontsize

#Run c++ program for results
os.system("g++ -std=c++11 -O3 %s -o proj && ./proj %s %d %g %g %g %d %d %d %s" % \
( programName, filename, L, T_init, T_final, h_T, mc_init, mc_final, h_mc, up ) )

#Extract results
data = np.loadtxt(filename + str(L), unpack=True)

#Unloading results
E_means = data[1]
M_means = data[-3]
MCcycles = np.log10(data[-1])

#Expectation value of energy as function of #MC-cycles
plt.figure(figsize=(10,8))
plt.title("L = %d, T = %g" % (L, T_final) )
plt.plot(MCcycles, E_means, label=r"$\langle E \rangle$")
plt.axhline(y=E_means[0], color='r', linestyle='--', label="Initial State Energy")
plt.xlabel(r"$\log_{10}\left(\# MC-cycles\right)$", fontsize=F)
plt.ylabel(r"$\langle E \rangle$", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.legend( fontsize=F)
plt.grid(True)

#Expectation value of magnetization as function of #MC-cycles
plt.figure(figsize=(10,8))
plt.title("L = %d, T = %g" % (L, T_final) )
plt.plot(MCcycles, M_means, label=r"$\langle \left| M \right| \rangle$")
plt.axhline(y=M_means[0], color='r', linestyle='--', label="Initial State Magnetization")
plt.xlabel(r"$\log_{10}\left(\# MC-cycles\right)$", fontsize=F)
plt.ylabel(r"$\langle \left| \mathcal{M} \right| \rangle$", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.legend( fontsize=F)
plt.grid(True)

#Number of accepted spin flips as function of MC-cycles
plt.figure(figsize=(10,8))
plt.title("L = %d, T = %g" % (L, T_final) )
plt.plot(MCcycles, data[-2], label=r"$\# N_A$")
plt.xlabel(r"$\log_{10}\left(\# MC-cycles\right)$", fontsize=F)
plt.ylabel(r"$\# N_A$", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.legend( fontsize=F)
plt.grid(True)

plt.show()
