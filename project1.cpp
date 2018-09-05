#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
using namespace std;
inline double f(double x){return 100.0*exp(-10.0*x);}


int main(int argc, char* argv[])
  {

    int n = atoi(argv[1]); //input number of grid points

    //Allocating memory for 6 vectors; a, b, c, b_tilde, d_tilde and v (solution)
    double* a       = new double[n-1];
    double* d       = new double[n];
    double* c       = new double[n-1];

    double* b       = new double[n];
    double* b_tilde = new double[n];
    double* d_tilde = new double[n];
    double* v       = new double[n];

    double h = (double)1/(n);
    double hh = h*h;

    for(int i = 0; i < n; i++)
      {
        a[i] = -1.0;
        d[i] = 2.0;
        c[i] = -1.0;
      }

      //Setting d_tilde equal to 2
      d_tilde[0] = d[0];
      b_tilde[0] = 0.0;     //Initial value of h^2 * f(x=0)=0


      /*Algorithm*/

      /*Forward substitution*/
    for(int i = 1; i <= n; i++)
      {
        d_tilde[i] = d[i] - (a[i-1]*c[i-1])/d_tilde[i-1];
        b_tilde[i] = f(i*h)*hh - (a[i-1]*b_tilde[i-1])/d_tilde[i-1];
      }

//Slette gir feilmelding delete[] a;

      v[n-1] = b_tilde[n-1]/d_tilde[n-1];

      /*Backward substitution*/

    for(int i=n-2 ; i >= 1; i--)
      {
        v[i] = (b_tilde[i] - c[i]*v[i+1])/d_tilde[i];
      }

    ofstream outfile;
    outfile.open("proj1aRun.dat");
    for(int i = 0; i <= n; i++)
      {
      outfile << i*h << " " << v[i] << endl;
      }


    return 0;
  }
