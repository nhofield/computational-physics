#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
using namespace std;
inline double f(double x){return 100.0*exp(-10.0*x);}
inline double analytical(double x){return 1.0-(1-exp(-10.0))*x - exp(-10.0*x);}

/*double norm(double*& array[], int sizeofarray)
  Calculates norm of a vector */ /*
  {
    double norm_computed = 0.0;
    for(int i=0 ; i<sizeofarray; i++)
      {
        norm_computed += array[i]*array[i]);
      }
    return sqrt(norm_computed);
  }
*/

/*
double subtract_arrays(double arr1[], double arr2[], int sizeofarrays)
  {
    double difference_array[sizeofarrays];
    for(int i=0; i < sizeofarrays; i++)
      {
        difference_array[i] = arr1[i]-arr2[i];
      }
    return difference_array;
  }
*/
int main(int argc, char* argv[])
  {

    int n = atoi(argv[1]); //input number of grid points

    //Allocating memory for 6 vectors; a, b, c, b_tilde, d_tilde and v (solution)
    double* a       = new double[n-1];
    double* d       = new double[n];
    double* c       = new double[n-1];

    double* u       = new double[n];  //Analytical solution
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

        //Analytical solution
        u[i] = analytical(h*i);
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

        //Analytical solution
      }

//Slette gir feilmelding delete[] a;

      v[n-1] = b_tilde[n-1]/d_tilde[n-1];

      /*Backward substitution*/

    for(int i=n-2 ; i >= 1; i--)
      {
        v[i] = (b_tilde[i] - c[i]*v[i+1])/d_tilde[i];
      }

      //Calc eps_relative
      double eps_relative;
      double temp_sum1 = 0.0;
      double temp_sum2 = 0.0;

      for(int i=0 ; i< n ; i++)
        {
          temp_sum1 += pow((u[i]-v[i]), 2);
          temp_sum2 += pow(u[i], 2);
        }
      eps_relative = sqrt(temp_sum1/temp_sum2);

      ofstream outfile;
      outfile.open("proj1aRun.dat");
      outfile << eps_relative << endl;
      outfile.close();
    return 0;
  }
