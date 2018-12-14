#include<iostream>
#include<algorithm>
#include<cstdlib>
#include<cmath>
#include<iterator>
using namespace std;
inline double N_rot(int N) {return (3.0/2)*N*(N-1);}

bool CheckIfContained(double* &Arr, double val, int N)
  {
    bool answer;
    int k = 0;
    double tol = 1.0e-14;

    for(int i = 0; i < N; i++)
      {
        if( fabs(val - Arr[i]) < tol)
          {
            answer = 2>1;
            break;
          }
        else
          {
            k += 1;
          }
      }
    if(k == N)
      {answer = 2<1;}
    else{}
    return answer;
  }


double* Sort_Array_Ascending(double* &ToBeSorted, int N)
  {
    double* SortedArr = new double[N];
    double* min_prev = new double[N];

    for(int i = 0; i < N; i++)
      {
        SortedArr[i] = ToBeSorted[i];

        for(int j = 0; j < N; j++)
          {
            if( ToBeSorted[j] < SortedArr[i] && CheckIfContained(ToBeSorted, ToBeSorted[j], N) == 0)
              {
                SortedArr[i] = ToBeSorted[j];
              }
            else{}
          }
      }
    return SortedArr;
  }

int main()
  {
    int n = 5;
    double* Messy_Array = new double[n];
    Messy_Array[0] = 3;
    Messy_Array[1] = 0;
    Messy_Array[2] = 100;
    Messy_Array[3] = -10;
    Messy_Array[4] = -2;
    double d = 10;
    double all = 0;
    double u = 100;
    double* Sorted_Array = Sort_Array_Ascending(Messy_Array, n);
    for(int i = 0; i < n; i++) cout << "Index = " << i << " Value = " << Sorted_Array[i] << endl;
    cout << CheckIfContained(Messy_Array, all, n) << endl;
    return 0;
  }
