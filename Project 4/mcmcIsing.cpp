#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include <chrono>
#include <ratio>
#include <ctime>

using namespace std;

double** MatrixAlloc(int);
void MatrixDeAlloc(double** &, int);
void MetropolisSampling(int, int, double, double* &, string);
void InitializeLattice(int, double** &, double &, double &, string);
void output(int, int, double, double* &);

// output file
ofstream outfile;

// inline function for periodic boundary conditions
inline int periodic(int i, int latticeLength, int add){return (i + add + latticeLength) % (latticeLength);}


int main(int argc, char* argv[])
  {
    string filename;  //String filename for storing output
    int L, mc_init, mc_final, h_mc;  //L: length of lattice [atoms]
    double T_init, T_final, h_T;  //h_T: step-size for temperature iteration
    string up;    //Specified 'u' for initializing lattice with all spins pointing up

    if ( argc > 1 )
      {
        filename =      argv[1];
        L        = atoi(argv[2]);
        T_init   = atof(argv[3]);
        T_final  = atof(argv[4]);
        h_T      = atof(argv[5]);
        mc_init  = pow(10, atoi(argv[6]) );
        mc_final = pow(10, atoi(argv[7]) );
        h_mc     = atoi(argv[8]);
        up       = argv[9];
      }
    // Declare new file name and add lattice size to file name
    outfile.open(filename.append(to_string(L)));
//    outfile << "    Temperature   |   <E>   |     C_V      |    <M>    |     X     |     |<M>|    |"  << endl;
    // Start Monte Carlo sampling by looping over T first
    int runNumber = 1;
    double* StoreValues = new double[6];
    for(double Temperature = T_init; Temperature <= T_final; Temperature += h_T)
      {
        for(int MCcycles = mc_init; MCcycles <= mc_final; MCcycles += h_mc )
          {
            cout << "Run " << runNumber << endl;

            for(int i = 0; i < 6; i++) StoreValues[i] = 0;
            // Start Monte Carlo computation
            auto start = chrono::high_resolution_clock::now();
            MetropolisSampling(L, MCcycles, Temperature, StoreValues, up);
            auto stop = chrono::high_resolution_clock::now();
            auto diff = stop-start;
            cout  << "Runtime of algorithm = " << chrono::duration <double,milli> (diff).count() << "ms" << endl;

            // Write current output to file
            output(L, MCcycles, Temperature, StoreValues);


            runNumber += 1;
          }
      }
    delete[] StoreValues;
    outfile.close(); // close output file

    return 0;
  }



void MetropolisSampling(int latticeLength, int MCcycles, double T, double* &StoreValues, string up)
  {
    // Initialize the seed and call the Mersenne algo
    std::random_device rd;
    std::mt19937_64 gen(rd());

    // Then set up the uniform distribution for x \in [0, 1]
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Allocate memory for lattice matrix
    double** Lattice = MatrixAlloc(latticeLength);

    // initialise energy and magnetization
    double Energy = 0.; double MagneticMoment = 0.;

    // initialize array for expectation values
    InitializeLattice(latticeLength, Lattice, Energy, MagneticMoment, up);
    // setup array for possible energy changes
    int E_max = 2*latticeLength*latticeLength;
    int vecLen = 2*E_max + 1;

    double* EnergyDifference = new double[vecLen];

    for(int dE =-E_max; dE <= E_max; dE += 4)
      {
        EnergyDifference[dE + E_max] = exp(-dE/T);
      }


    int Accepted = 0; //To store number of accepted flips

    for (int cycle = 1; cycle <= MCcycles; cycle++)
      {
        // The sweep over the lattice, looping over all spin sites
        for(int x = 0; x < latticeLength; x++)
          {

            for(int y = 0; y < latticeLength; y++)
              {
                int ix = (int) (distribution(gen)*(double)latticeLength);
                int iy = (int) (distribution(gen)*(double)latticeLength);
                int deltaE = 2*Lattice[ix][iy]*
                ( Lattice[ix][periodic(iy, latticeLength,-1)] + Lattice[periodic(ix, latticeLength,-1)][iy]
                + Lattice[ix][periodic(iy, latticeLength, 1)] + Lattice[periodic(ix, latticeLength, 1)][iy] );

                if ( deltaE <= 0 )
                  {
                    Lattice[ix][iy] *= -1.0; // flip one spin and accept new spin config
                    MagneticMoment += (double) 2*Lattice[ix][iy];
                    Energy += (double) deltaE;
                    Accepted += 1;
                  }


               else
                {
                  if ( distribution(gen) <= EnergyDifference[deltaE + E_max] )
                    {
                      Lattice[ix][iy] *= -1.0; // flip one spin and accept new spin config
                      MagneticMoment += (double) 2*Lattice[ix][iy];
                      Energy += (double) deltaE;
                      Accepted += 1;
                    }
                  else
                    {
                      //Do nothing
                    }
                }

              }
          }

        // update expectation values for local node
        StoreValues[0] += Energy; StoreValues[1] += Energy*Energy;
        StoreValues[2] += MagneticMoment;
        StoreValues[3] += MagneticMoment*MagneticMoment;
        StoreValues[4] += fabs(MagneticMoment);
        StoreValues[5] = Accepted;
      }
    MatrixDeAlloc(Lattice, latticeLength);    //Deallocate memory
    delete[] EnergyDifference;
  } // end of Metropolis sampling over spins

void InitializeLattice(int latticeLength, double** &Lattice, double &Energy, double &MagneticMoment, string up)
// setup spin matrix and initial magnetization
  {
    if( up == "u")
    //Set all spin orientations up (+1)
      {
        for(int x = 0; x < latticeLength; x++)
          {
            for (int y= 0; y < latticeLength; y++)
              {
                Lattice[x][y] = 1.0; // spin orientation for the ground state
                MagneticMoment += (double) Lattice[x][y];
              }
          }
      }
    else
    //Set all spin orientations randomly
      {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        int SpinPossibilities[2] = {-1, 1};

        for(int x = 0; x < latticeLength; x++)
          {
            for (int y= 0; y < latticeLength; y++)
              {
                int randSpin = (int) (distribution(gen)*(double)2);
                Lattice[x][y] = SpinPossibilities[randSpin]; // spin orientation for the ground state
                MagneticMoment += (double) Lattice[x][y];
              }
          }
      }
      // setup initial energy
    for(int x = 0; x < latticeLength; x++)
      {
        for (int y = 0; y < latticeLength; y++)
          {
            Energy -= (double) Lattice[x][y]*(Lattice[periodic(x, latticeLength,-1)][y] +
                                        Lattice[x][periodic(y, latticeLength,-1)] );
          }
      }

  }// end function InitializeLattice

void output(int latticeLength, int MCcycles, double T, double* &StoreValues)
  {

    double norm = 1.0/((double) (MCcycles)); // divided by number of cycles
    double E_ExpectationValues = StoreValues[0]*norm;
    double E2_ExpectationValues = StoreValues[1]*norm;
    double M_ExpectationValues = StoreValues[2]*norm;
    double M2_ExpectationValues = StoreValues[3]*norm;
    double Mabs_ExpectationValues = StoreValues[4]*norm;
    // all expectation values are per spin, divide by 1/NSpins/NSpins
    double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/latticeLength/latticeLength;
    double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/latticeLength/latticeLength;
    outfile << setiosflags(ios::showpoint | ios::uppercase);
    outfile << setw(15) << setprecision(8) << T;
    outfile << setw(15) << setprecision(8) << E_ExpectationValues/latticeLength/latticeLength;
    outfile << setw(15) << setprecision(8) << Evariance/T/T;
    outfile << setw(15) << setprecision(8) << M_ExpectationValues/latticeLength/latticeLength;
    outfile << setw(15) << setprecision(8) << Mvariance/T;
    outfile << setw(15) << setprecision(8) << Mabs_ExpectationValues/latticeLength/latticeLength;
    outfile << setw(15) << setprecision(8) << StoreValues[5];
    outfile << setw(15) << setprecision(8) << MCcycles << endl;
  } // end output function


double** MatrixAlloc(int n)
/*Function which returns dynamically allocates memory for a square nxn matrix*/
  {
    double** Matrix = new double*[n];
    for(int i = 0; i < n; i++)
      {
        Matrix[i] = new double[n];
      }
    return Matrix;
  }
void MatrixDeAlloc(double** &Matrix, int size)
/*Void for deallocating dynamic matrix*/
  {
    for(int i = 0; i < size; i++) delete[] Matrix[i];
    delete[] Matrix;
  }
