#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <chrono>
#include <ratio>
#include <ctime>
#include <algorithm>

using namespace std;
inline double f(double x){return 100.0*exp(-10.0*x);}
inline double analytical(double x){return 1.0-(1-exp(-10.0))*x - exp(-10.0*x);}

int main(int argc, char* argv[])
  {

    int n = atoi(argv[1]); //input number of grid points

    //Allocating memory for 6 vectors; a, b, c, b_tilde, d_tilde and v (solution)
    double* u       = new double[n+2];  //Analytical solution
    double* b_tilde = new double[n];
    double* d_tilde = new double[n];
    double* v       = new double[n+2];
    double* b       = new double[n];

    double h = (double)1/(n+1);
    double hh = h*h;

    //Computing analytical solution
    for(int i = 1; i <= n; i++)
      {
        u[i] = analytical(i*h);
      }

    for(int i = 0; i < n; i++)
      {
        b[i] = hh*f(h*(i+1));
      }


    //Boundary conditions
    u[0] = u[n+1] = 0.0;
    v[0] = v[n+1] = 0.0;

    //b_tilde and d_tilde initial conditions
    b_tilde[0] = b[0];
    d_tilde[0] = 2.0;

    /*START Algorithm*/

    auto start = chrono::high_resolution_clock::now();

    /*Forward substitution*/
    for(int i = 1; i < n; i++)
      {
        d_tilde[i] = (double) (i+2)/(i+1);
        b_tilde[i] = b[i] + (double) i/(i+1) * b_tilde[i-1];
      }

    //Compute v_n
    v[n] = b_tilde[n-1]/d_tilde[n-1];

      /*Backward substitution*/
    for(int i = n-1 ; i >= 1; i--)
      {
        v[i] = (double) i/(i+1) * (b_tilde[i-1] + v[i+1]);
      }

    /*STOP Algorithm*/

    auto stop = chrono::high_resolution_clock::now();
    auto diff = stop-start;
    cout  << " Runtime of algorithm = " << chrono::duration <double,milli> (diff).count() << "ms" << endl;

    //Not used anymore
    delete []b_tilde;
    delete []d_tilde;
    delete []b;

    //Calc eps_relative

    //Writing result to file
    ofstream outfile;
    outfile.open("proj1data.dat");

    if(*argv[2] == 'Y')
      {
        for(int i = 1; i <= n; i++) outfile << i*h << " " << u[i] << " " << v[i] << endl;
        delete []u;
        delete []v;
      }
    else
      {
        //Alocating memory for storing relative error
        double* eps_relative = new double[n];
        for(int i = 0; i <= n; i++) eps_relative[i] = log10( abs( ( v[i+1] - u[i+1] )/u[i+1] ) );

        delete []u;
        delete []v;
        //Fetching maximum value in array, then deallocating memory
        outfile << *max_element(eps_relative, eps_relative + n) << endl;
        delete []eps_relative;
      }
    outfile.close();

    return 0;
  }
