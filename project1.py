import numpy as np
import os
import time
import matplotlib.pyplot as plt
"""
Author: Noah Oldfield
Program: Has two query inputs:
first specifies whether to run general, special or LU-decomposition algorithm.
second specifies whether to plot approx to analytical solutions in same plot
or plot the relative error.
Depends on the following programs:
project1.cpp, project1_imp.cpp, project1_LU.cpp
"""
N = 2
eps_relative = np.zeros(N)
number_N = np.zeros(N)

#First query
Q1 = raw_input("Type G for general algorithm, S for the special or LU for LU-decomposition: ")
print "---------------------------"
if Q1 == 'G':
    program = "project1.cpp"    #General algorithm
    dueToArmadillo = ""
elif Q1 == 'S':
    program = "project1_imp.cpp"    #Special algorithm
    dueToArmadillo = ""
elif Q1 == 'LU':
    program = "project1_LU.cpp"
    dueToArmadillo = "-O1 -larmadillo"

filename = "proj1data.dat"  #The file that the c++ program generates
print program
#Second query
Q2 = raw_input("""
Type Y to:
Compare numerical solution with analytical solution in plot.
-------------------------
Type ENTER/anything to:
Create log-log plot of relative error as a function of grid points.
-------------------------
""")
F = 14  #Fontsize in plot

for i in range(N):

    n = 10**(i+1)
    number_N[i] = n

    if Q2 == "Y":
        run = os.system("g++ -std=c++11 %s %s -o project1 && ./project1 %d Y && rm project1" % (program, dueToArmadillo, n))
        x, u, v = np.loadtxt(filename, unpack=True)
    else:
        run = os.system("g++ -std=c++11 %s %s -o project1 && ./project1 %d N && rm project1" % (program, dueToArmadillo, n))
        eps_relative[i] = np.loadtxt(filename, unpack=True)

plt.figure(figsize=(10,8))

if Q2 == "Y":
    plt.title("Comparison of analytical and numerical solutions.")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(x, u, label="Analytical solution u(x)")
    plt.plot(x, v, label="Numerical approximation v(x)")

else:
    plt.title("log-log plot of relative error as a function of number of grid points.")
    plt.xlabel(r"$log_{10}(n)$", fontsize=F)
    plt.ylabel(r"$log_{10}\epsilon_{relative}$", fontsize=F)
    plt.plot(np.log10(number_N), eps_relative, label=r"$log_{10}\epsilon_{relative}(n)$")

plt.grid(True)
plt.legend()
plt.show()
