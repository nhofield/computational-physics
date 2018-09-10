#include <iostream>
#include "armadillo"
#include <fstream>
#include <chrono>
#include <ratio>
#include <ctime>
#include <algorithm>
#include <cstdlib>

using namespace std;
using namespace arma;
inline double f(double x){return 100.0*exp(-10.0*x);}
inline double analytical(double x){return 1.0-(1-exp(-10.0))*x - exp(-10.0*x);}

int main(int argc, char* argv[])
  {
    int n = atoi(argv[1]); //input number of grid points

    //Defining matrix A and vectors b_tilde and u
    mat A(n,n);
    vec b_tilde(n, 1);
    vec u(n,1);

    //Populating A
    for(int i = 0; i < n; i++) A(i,i) = 2.0;
    for(int i = 0; i < n-1; i++) A(i+1, i) = A(i, i+1) = -1.0;

    double h = (double) 1/(n+1);
    double hh = h*h;

    //Analytical solution
    for(int i = 0; i < n; i++) u(i) = analytical(h*(i+1));
    //Filling b_tilde
    for(int i = 0; i < n; i++) b_tilde(i) = hh*f(h*(i+1));


    auto start = chrono::high_resolution_clock::now();
    /*START algorithm*/
    vec v = solve(A, b_tilde, solve_opts::fast);
    /*STOP Algorithm */
    auto stop = chrono::high_resolution_clock::now();
    auto diff = stop-start;
    cout  << " Runtime of algorithm = " << chrono::duration <double,milli> (diff).count() << "ms" << endl;

    A.reset();

    ofstream outfile;
    outfile.open("proj1data.dat");

    if(*argv[2] == 'Y')
      {
        for(int i = 0; i < n; i++) outfile << (i+1)*h << " " << u(i) << " " << v(i) << endl;
      }
    else
      {
        //Alocating memory for storing relative error
        double* eps_relative = new double[n];
        for(int i = 0; i < n; i++) eps_relative[i] = log10( abs( ( v(i) - u(i) )/u(i) ) );

        //Deallocating memory
        v.reset();
        u.reset();

        //Fetching maximum value in array, then deallocating memory
        outfile << *max_element(eps_relative, eps_relative + n) << endl;
        delete []eps_relative;
      }
    outfile.close();

    return 0;
  }
