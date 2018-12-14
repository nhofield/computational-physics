import numpy as np
import os
import matplotlib.pyplot as plt

"""Computes the analytical expressions for the 2x2 lattice for T=1:
  Z    : partition function
 <E>   : Expectation value of energy
<|M|>  : Expectation value of absolute magnetization
 C_V   : Specific heat
  X    : Magnetic susceptibility
"""
L = 2.
N = L**2
Z = 2.*( np.e**(8) + np.e**(-8) + 6)            #Z

E_mean_tot = 16.*( np.e**(-8) - np.e**8)/Z           #<E>
M_mean_tot = 8.*(np.e**8 + 2)/Z                      #<|M|>
E2_mean = 128.*( np.e**8 + np.e**(-8) )/Z
M2_mean = 32.*(np.e**8 + 1)/Z
C_V = (E2_mean - E_mean_tot**2)/N
X = (M2_mean - M_mean_tot**2)/N
E_mean = E_mean_tot/N                                #Per spin
M_mean = M_mean_tot/N                                #Per spin

print "Units per spin"
print "     <E>     ", "    <|M|>   ", "      C_V      ", "       X        "
print E_mean, M_mean, C_V, X

programName = "mcmcIsing.cpp"
filename = "results"

MCcycles = np.array([10**i for i in range(7)])
MClen = len(MCcycles)
E_mean_diff = np.zeros(MClen)
M_mean_diff = np.zeros_like(E_mean_diff)
C_V_diff = np.zeros_like(E_mean_diff)
X_diff = np.zeros_like(E_mean_diff)

F = 12

plt.figure()
for i in range(MClen):

    os.system("g++ -std=c++11 %s -o proj && ./proj %s %d %d 1 1 1 && rm proj" % (programName, filename, L, MCcycles[i]))
    data = np.loadtxt(filename + "2", unpack=True)
    E_mean_diff[i] = abs(E_mean - data[1])
    M_mean_diff[i] = abs(M_mean - data[-1])
    C_V_diff[i] = abs(C_V - data[2])
    X_diff[i] = abs(X - data[4])

MC_plot = np.log10(MCcycles)
plt.plot(MC_plot, np.log10(E_mean_diff), label=r"$log10\epsilon_{<E>}$")
plt.plot(MC_plot, np.log10(M_mean_diff), label=r"$log10\epsilon_{<|M|>}$")
plt.plot(MC_plot, np.log10(C_V_diff), label=r"$log10\epsilon_{C_V}$")
plt.plot(MC_plot, np.log10(X_diff), label=r"$log10\epsilon_{\chi}$")
plt.xlabel("log10(#MC-cycles)", fontsize=F)
plt.ylabel(r"$log10\epsilon$", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.legend( fontsize=F)
plt.grid(True)
plt.show()
