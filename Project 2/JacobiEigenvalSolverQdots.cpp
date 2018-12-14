#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <ratio>
#include <ctime>
#include <fstream>
#include "armadillo"

using namespace std;
using namespace arma;
#include "JacobiEigenvalSolver.h"
inline double N_rot(int N) {return (double)(3.0/2)*N*(N-1);}

const double pi = 3.14159265358979323846264338327950288419716939937;



int main()
  {
    /*Test calls*/
    test_Max_Off_Diag_kl();
    test_Max_Off_Diag_Element();
    test_JacobiRotation();
    test_JacobiEigenValueAlgorithm();
    test_dotproduct();
    test_FrobeniusNorm();
    cout << "------------------\n| Tests complete |\n------------------\n" << endl;
    cout << "Main-Program\n---------" << endl;

    //Size of matrix (nxn)
    const int N = 100;
    double rho_N = 10.0;
    double h = rho_N/N;
    double hh = h*h;
    double w_r = 5.0;
    double w_r_pow2 = w_r*w_r;

    //Tolerance for jacobi method when checking frobenius norm of off(A)
    double eps = 1.0e-16;

    //Defining A, and matrix V for storing eigenvectors
    double** A = MatrixAlloc(N);
    double** V = MatrixAlloc(N);

    //Diagonal elements of A
    double* d = new double [N];
    double* b = new double[N-1];

    for(int i = 0; i < N; i++)
      {
        d[i] = 2.0/hh + w_r_pow2*((i+1)*h)*((i+1)*h) + 1.0/((i+1)*h);  //Adding term to diagonal
      }
    for(int i = 0; i < N-1; i++) b[i] = -1.0/hh;

    //Creates the toeplitz matrix A
    CreateToeplitzMatrix(A, N, d, b);
    //Sets eigenvector matrix to identity matrix
    Identity(V, N);

    /*START Algorithm*/
    auto start = chrono::high_resolution_clock::now();
    double* EigenValues = JacobiEigenValueAlgorithm(A, V, eps, N);
    auto stop = chrono::high_resolution_clock::now();
    /*STOP Algorithm*/

    auto diff = stop-start;
    cout  << "Runtime of algorithm = " << chrono::duration <double,milli> (diff).count() << "ms" << endl;

    MatrixDeAlloc(A, N);
    MatrixDeAlloc(V, N);
    double* Analytical_Eig = Analytical_Eigenvalues_Schr(N);
    /*
    CompareEigenvalueComputations(Analytical_Eig, EigenValues, N);
    cout << "Eigenvalues\n----------" << endl;
    */
    PrintArray(EigenValues, N);

    /*
    char Q;
    cout << "Write Eigenvalues to file? (Y for yes or N to avoid): " << endl;
    cin >> Q;
    if(Q == 'Y')
      {*/
        ofstream outfile;
        outfile.open ("eigenvaluesQdots.dat");
        for(int i = 0; i < N; i++)
          {
            outfile << setprecision(8) << EigenValues[i] << " " << Analytical_Eig[i] << endl;

          }
        outfile.close();
        outfile.open("eigenvectorQdots.dat");
        for(int i = 0; i < N; i++)
          {
            outfile << setprecision(8) << V[0][i] << " " << (i+1)*h << endl;
          }
        outfile.close();
/*      }
    else{}
*/

    return 0;
  }

double* Analytical_Eigenvalues_Schr(int N)
    /*The function for analytical eigenvalues of buckling-beam problem*/
    {
      double* Analytical_Eig = new double[N];
      for(int i = 0; i < N; i++)
        {
          Analytical_Eig[i] = 2.0*i + 3.0/2.0 ;
        }
      return Analytical_Eig;
    }

double Max_Off_Diag_Element(double** &A, int n)
    /*Finds absolute value of maximum non diagonal element of given matrix A*/
    {
      double A_max_offdiag = 0.0;

      for(int i = 0; i < n; i++)
        {
          for(int j = i+1; j < n; j++)
            {
              if( fabs(A[i][j]) > A_max_offdiag )
                {
                  A_max_offdiag = fabs(A[i][j]);
                }
              else
                {

                }
            }
        }

      return A_max_offdiag;
    }
void test_Max_Off_Diag_Element()
      /**/
    {
      int N = 4;
      double** A_test = MatrixAlloc(N);

      for(int i = 0; i < N; i++) A_test[i][i] = A_test[0][i] = A_test[i][0] = i;
      A_test[2][1] = A_test[1][2] = A_test[3][1] = A_test[1][3] = -10;
      A_test[3][2] = A_test[2][3] = 100;

      double ActualMaxElement_A_nondiag = 100;
      double tol = 1.0e-16;
      double val;
      val = Max_Off_Diag_Element(A_test, N);
      if( fabs( val - ActualMaxElement_A_nondiag) < tol )
        {

        }
      else
        {
          cout << "Error in Max_Off_Diag_Element, Actual = " << ActualMaxElement_A_nondiag << endl;
          cout << "Function gives " << val << endl;
        }
      MatrixDeAlloc(A_test, N);
    }

double dotproduct(double* &Arr1, double* &Arr2, int n)
      {
        double dotprod_result;
        dotprod_result = 0.0;

        for(int i = 0; i < n; i++)
          {
            dotprod_result += Arr1[i]*Arr2[i];
          }
        return dotprod_result;
      }
void test_dotproduct()
      {
        int N = 4;
        double tol = 1e-16;
        double* a = new double[N];
        double* b = new double[N];
        for(int i = 0; i < N; i++)
          {
            a[i] = i;
            b[i] = i-1;
          }
        double answer = 8;
        double test_result = dotproduct(a, b, N);
        if( fabs(answer - test_result) < tol )
          {

          }
        else
        {
          cout << "Error in dotproduct function, Actual = " << answer << endl;
          cout << "Function returns = " << test_result << endl;
        }
      }

void Max_Off_Diag_kl(double** &A, int* k, int* l, int n)
      /*Finds indices of absolute value of maximum element of given matrix*/
    {
      *k = 0;
      *l = 1;
      for(int i = 0; i < n; i++)
      {
        for(int j = i+1; j < n; j++)
        {
          if( fabs(A[i][j]) > fabs(A[*k][*l]) )
          {
            *k = i;
            *l = j;
          }
          else
          {

          }//end if
        }//end for
      }//end for
    }
void test_Max_Off_Diag_kl()
    {
      int N = 3;
      double** A_test = MatrixAlloc(N);
      for(int i = 0; i < N; i++) A_test[0][i] = A_test[i][0] = i+1;
      A_test[1][1] = 4.0;
      A_test[2][2] = 8.0;
      A_test[1][2] = A_test[2][1] = 5.0;
      int ActualMax_k = 1;
      int ActualMax_l = 2;
      int p, q;
      double tol = 1.0e-16;
      Max_Off_Diag_kl(A_test, &p, &q, N);
      if( fabs(p - ActualMax_k) < tol && fabs(q - ActualMax_l) < tol)
      {

      }
      else
      {
        cout << "Error in Max_Off_Diag_kl function" << "Actual = (" << ActualMax_k << "," << ActualMax_l << ")" << endl;
        cout << "Function returns (" << p << "," << q << ")" << endl;
      }
      MatrixDeAlloc(A_test, N);
    }

double* JacobiRotation(double** &A, int k, int l)
      /*Performs a jacobi rotation of given matrix A*/
    {
      double tau, t_min;
      int c_s_len = 2;
      double* c_s = new double[c_s_len];

      if( fabs( A[k][l] ) != 0.0)
      {
        tau = (A[l][l] - A[k][k])/(2.0*A[k][l]);
        if( tau >= 0.0 )
        {
          t_min =  - tau + sqrt(1.0 + tau*tau);
        }
        else
        {
          t_min = - ( tau + sqrt(1.0 + tau*tau) );
        }
        c_s[0] = 1.0/sqrt(1.0 + t_min*t_min);
        c_s[1] = t_min*c_s[0];
      }
      else
      {
        c_s[0] = 1.0;
        c_s[1] = 0.0;
      }
      return c_s;
    }
void test_JacobiRotation()
    {
      int N = 2;
      double tol = 1.0e-16;
      double** A_test = MatrixAlloc(N);
      A_test[0][0] = A_test[1][1] = 1.0;
      A_test[0][1] = A_test[1][0] = 2.0;
      int k, l;
      k = 0;
      l = 1;
      double* result = JacobiRotation(A_test, k, l);
      double actual_c, actual_s;
      actual_c = actual_s = 1.0/sqrt(2);
      if( fabs(result[0] - actual_c) < tol && fabs(result[1] - actual_s) < tol)
      {

      }
      else
      {
        cout << "Error in Jacobi_rotation, actual c, s = " << actual_c <<"," << actual_s << endl;
        cout << "Function returns " << result[0] << "," << result[1] << endl;
      }
      MatrixDeAlloc(A_test, N);
    }

void Update_A_V(double** &A, double** &V, int n, int k, int l, double* &c_s)
        /*Performs orthogonal transformation with jacobi rotation matrix*/
      {
        double a_kk, a_ll, a_ik, a_il, c, s;

        c = c_s[0];
        s = c_s[1];

        a_kk = A[k][k];
        a_ll = A[l][l];

        A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
        A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
        A[k][l] = A[l][k] = 0.0;

        for(int i = 0; i < n ; i++)
          {
            if( i != k && i != l)
              {
                a_ik = A[i][k];
                a_il = A[i][l];
                A[i][k] = c*a_ik - s*a_il;
                A[k][i] = A[i][k];
                A[i][l] = c*a_il + s*a_ik;
                A[l][i] = A[i][l];
              }
            else
              {

              }
            double V_ik, V_il;
            V_ik = V[i][k];
            V_il = V[i][l];
            V[i][k] = c*V_ik - s*V_il;
            V[i][l] = c*V_il + s*V_ik;
          }
      }

bool Test_Orthogonality(double** &V, double* NowFrobeniusNorm, int N)
      {
        double WhatNowFrobeniusNorm = FrobeniusNorm(V, N);
        double tol = 1.0e-11;

        return fabs(*NowFrobeniusNorm-WhatNowFrobeniusNorm) < tol;
      }

double* JacobiEigenValueAlgorithm(double** &A, double** &V, double eps, int n)
          /*Returns array of eigenvalues of symmetric matrix A*/
        {
          double* EigenValues = new double[n];
          //Compute initial frobenius norm of V
          double FrobeniusNorm_V;
          FrobeniusNorm_V = FrobeniusNorm(V, n);

          int k, l, c, s;
          int counter = 0;
          //Skaff nåværende k og l slik max|a_ij|
          while(Max_Off_Diag_Element(A, n) > eps)
            {

              if( Test_Orthogonality(V, &FrobeniusNorm_V, n) != 1 )
                {
                  cout << "Orthogonality not conserved in the " << counter << "rotation. Loop not broken, eigenvalues possibly not correct" << endl;
                  break;
                }
              else //Orthogonality conserved, keep going
              {
                Max_Off_Diag_kl(A, &k, &l, n);  //Changes to current k,l <-> max|aij|
                double* c_s = JacobiRotation(A, k, l);
                Update_A_V(A, V, n, k, l, c_s);
                counter += 1;
              }

            }

          //Retrieve eigenvalues from diagonal
          for(int i = 0; i < n; i++) EigenValues[i] = A[i][i];
          cout << counter << " Rotations" << endl;
          cout << "Expected Number of rotations " << N_rot(n) << endl;
          return EigenValues;
        }
void test_JacobiEigenValueAlgorithm()
          {
            int n = 3;
            double** A_test = MatrixAlloc(n);
            double** V_test = MatrixAlloc(n);
            Identity(V_test, n);
            A_test[0][0] = A_test[0][1] = A_test[1][0] = A_test[1][2] = A_test[2][1] = 1.0;
            A_test[0][2] = A_test[2][0] = -2.0;
            A_test[1][1] = 2.0;
            A_test[2][2] = -1.0;
            double Actual1 = 2.0;
            double Actual2 = sqrt(7);
            double Actual3 = -sqrt(7);
            double eps = 1.0e-9;  //for jacobi method
            double tol = 1.0e-14;  //for equality test
            PrintMatrix(A_test, n);
            double* result = JacobiEigenValueAlgorithm(A_test, V_test, eps, n);
            if(fabs(result[0]-Actual1) < tol && fabs(result[1] - Actual2) < tol && fabs(result[2]-Actual3)<tol)
              {

              }
            else
              {
                cout << "Error in JacobiEigenValueAlgorithm, correct = " << Actual1 << "," << Actual2 << "," << Actual3 << ")" << endl;
                cout << "Returns (" << result[0] << "," << result[1] << "," << result[2] << ")" << endl;
              }
          MatrixDeAlloc(A_test, n);
          MatrixDeAlloc(V_test, n);
          }

void PrintMatrix(double** &Matr, int n)
    {
      for(int i = 0; i < n ; i++)
        {
          for(int j = 0; j < n; j++) cout << setw(7) << Matr[i][j] << " ";
          cout << endl;
        }
    }
void PrintArray(double* &Arr, int n)
  {

    cout << "Array=" << endl;
    for(int i = 0; i < n; i++)
      {
        cout << Arr[i] << endl;
      }
  }

double FrobeniusNorm(double** &Matr, int N)
    {
      double frobeniusNorm_squared = 0.0;
      double frobeniusNorm = 0.0;

      for(int i = 0; i < N; i++)
        {
          for(int j = 0; j < N; j++)
            {
              frobeniusNorm_squared += Matr[i][j]*Matr[i][j];
            }
        }
      frobeniusNorm = sqrt(frobeniusNorm_squared);
      return frobeniusNorm;
    }
void test_FrobeniusNorm()
    /*Two test cases, frobenius norm of identity matrix of size n should be sqrt(n)
    and another matrix*/
    {
      int N = 3;
      double tol = 1.0e-16;
      double** V_test1 = MatrixAlloc(N);
      double** V_test2 = MatrixAlloc(N);

      for(int i = 0; i < N; i++)
        {
          V_test1[0][i] = 1.0;
          V_test1[1][i] = 2.0;
          V_test1[2][i] = 3.0;
        }
      Identity(V_test2, N);

      double ActualFrobeniusNorm1 = sqrt(42);
      double ActualFrobeniusNorm2 = sqrt(N);
      double ComputedFrobeniusNorm1 = FrobeniusNorm(V_test1, N);
      double ComputedFrobeniusNorm2 = FrobeniusNorm(V_test2, N);

      if(fabs(ActualFrobeniusNorm1 - ComputedFrobeniusNorm1) < tol && fabs(ActualFrobeniusNorm2 - ComputedFrobeniusNorm2) < tol)
        {

        }
      else
        {
          cout << "Error in FrobeniusNorm function, Actual1 = " << ActualFrobeniusNorm1 << " and Actual2 = " << ActualFrobeniusNorm2 << endl;
          cout << "Function returns " << ComputedFrobeniusNorm1 << " and " << ComputedFrobeniusNorm2 << endl;
        }
    }

double* Analytical_Eigenvalues_Buck(int N, double d, double a, double hh)
  /*The function for analytical eigenvalues of buckling-beam problem*/
  {
    double* Analytical_Eig = new double[N];
    for(int i = 0; i < N; i++)
      {
        Analytical_Eig[i] = (2.0/hh)*(1.0 - cos( (double)(pi/(N+1))*(i+1.0) ));
      }
    return Analytical_Eig;
  }
void CreateToeplitzMatrix(double** &A, int n, double* &d, double* &b)
  /*Creates Toeplitz matrix, d = diagonal elements, b = diagonal elements just
  above main diagonal and below main diagonal*/
  {
    for(int i = 0; i < n; i++)
      {
        for(int j = 0; j < n; j++)
          {
            A[i][j] = 0.0;
            if(i == j)
              {
                A[i][j] = d[i];
              }
            else{}
          }
      }
    for(int i = 0; i < n-1; i++) A[i][i+1] = A[i+1][i] = b[i];
  }
void Identity(double** &Matr, int n)
    /*Creates Identity matrix of size n*/
    {
      for(int i = 0; i < n; i++)
        {
          for(int j = 0; j < n; j++)
            {
              Matr[i][j] = 0.0;
              if(i==j)
                {
                  Matr[i][j] = 1.0;
                }
              else{}
            }
        }
    }
void CompareEigenvalueComputations(double* &Analytical, double* &JacobiComputed, int n)
    /*Takes two arrays of eigenvalues, one computed numerically and the other analytically.
      Since the eigenvalues in the arrays may not be in the same order, it tests all
      elements against one another.*/
    {
      int g = 0; // Counter
      double tol = 1.0e-8;
      for(int i = 0; i < n; i++)
        {
          for(int j = 0; j < n; j++)
            {
              if( fabs(Analytical[i] - JacobiComputed[j]) < tol )
                {
                  g += 1;
                }
              else{}

            }
        }
      if(fabs(g - n) > tol)
        {
          cout << "Error CompareEigenvalueComputations; Eigenvalues not equal or tolerance " << tol << " too low" << endl;
        }
      else
        {
          cout << "Success, eigenvalues are equal to analytical eigenvalues within tolerance " << tol << endl;
        }
    }

double** MatrixAlloc(int n)
    /*Function which returns dynamically allocates memory for a mxn matrix*/
      {
        double** Matrix = new double*[n];
        for(int i = 0; i < n; i++)
          {
            Matrix[i] = new double[n];
          }
        return Matrix;
      }
void MatrixDeAlloc(double** &Matrix1, int size)
    /*Void for deallocating dynamic matrix*/
      {
        for(int i = 0; i < size; i++) delete Matrix1[i];
        delete Matrix1;
      }
