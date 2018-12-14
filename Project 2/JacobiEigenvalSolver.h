#ifndef JACOBIEIGENVALSOLVER_H
#define JACOBIEIGENVALSOLVER_H

double    Max_Off_Diag_Element(double** &A, int n);
void test_Max_Off_Diag_Element();

void      Max_Off_Diag_kl(double** &A, int* k, int* l, int n);
void test_Max_Off_Diag_kl();

double* JacobiRotation(double** &A, int k, int l);
void test_JacobiRotation();

void Update_A_V(double** &A, int n, int k, int l, double* &c_s);

double* JacobiEigenValueAlgorithm(double** &A, double** &V, double eps, int n);
void test_JacobiEigenValueAlgorithm();

void PrintMatrix(double** &Matr, int n);
void PrintArray(double* &Arr, int n);

double* Analytical_Eigenvalues_Schr(int N);

void CreateToeplitzMatrix(double** &A, int n, double* &d, double* &b);
void Identity(double** &Matr, int n);

double dotproduct(double* &Arr1, double* &Arr2, int n);
void test_dotproduct();

double FrobeniusNorm(double** &Matr, int N);
void test_FrobeniusNorm();

double** MatrixAlloc(int n);
void MatrixDeAlloc(double** &Matrix1, int size);

double* Analytical_Eigenvalues_Buck(int N, double d, double a, double hh);
void CompareEigenvalueComputations(double* &Analytical, double* &JacobiComputed, int n);

#endif
