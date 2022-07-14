//==============================================================================
// Start of the file "hybrid_fields.cc"
//
// Version 3b.
//
// Written by Vladimir Florinski
//==============================================================================


#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include "geo_memory.hh"
#include "geo_coord.hh"
#include "block_thomas.hh"
#include "hybrid_fields.hh"

using namespace std;


//==============================================================================
// The fields_t class public methods
//==============================================================================


//------------------------------------------------------------------------------
// Allocate memory for fields arrays
//------------------------------------------------------------------------------
// Input  | cells    | number of grid cells
// Input  | length   | size of the box
// Input  | fratio   | plasma to cyclotron frequency ratio
//------------------------------------------------------------------------------

void fields_t::Activate(int cells, double length, double fratio)
{

// Can only be activated once
   if(isfld) return;

// Key parameters
   Imax = cells + 2; // 2 ghost cells
   xmax = length;
   dx = xmax / (Imax - 2);
   wpiwci = fratio;

// Memory for the EM fields and other bulk quantities
   B = Create2D<double>(Imax, 3);
   E = Create2D<double>(Imax, 3);
   j = Create2D<double>(Imax, 3);

   n = new double[Imax + 1];
   u = Create2D<double>(Imax, 3);
   p = new double[Imax + 1];

   isfld = true;
};


//------------------------------------------------------------------------------
// Class constructor with memory allocation
//------------------------------------------------------------------------------
// Input  | cells    | number of grid cells
// Input  | length   | size of the box
// Input  | fratio   | plasma to cyclotron frequency ratio
//------------------------------------------------------------------------------

fields_t::fields_t(int cells, double length, double fratio)
{
   isfld = false;
   Activate(cells, length, fratio);
};


//------------------------------------------------------------------------------
// Apply the initial condition for the EM field
//------------------------------------------------------------------------------
// Input  | type     | type of initial field (future work)
// Input  | B0[4]    | constant background magnetic field vector
// Input  | epre     | initial electron pressure
//------------------------------------------------------------------------------

void fields_t::Init(int type, double *B0, double epre)
{

   int i;

   if(!isfld) return;
   pe0 = epre;

// Generate the desired EMF
   for(i = 1; i <= Imax; i++) {
      switch(type) {

// Uniform magnetostatic
      case 1:
         memcpy(&B[i][1], &B0[1], 3 * sizeof(double));
         memset(&E[i][1], 0, 3 * sizeof(double));
         break;

// Zero field
      default:
         memset(&B[i][1], 0, 3 * sizeof(double));
         memset(&E[i][1], 0, 3 * sizeof(double));
         break;
      };
   };
};


//------------------------------------------------------------------------------
// Return the interpolated EM fields
//------------------------------------------------------------------------------
// Input  | x        | x coordinate
// Output | Bint[4]  | magnetic field vector
// Output | Eint[4]  | electric field vector
//------------------------------------------------------------------------------

void fields_t::GetEMF(double x, double *Bint, double *Eint)
{
   int xyz, i[3];
   double f[3], res = 0.0;

   if(!isfld) return;

// If the flag is set, we add the term 𝜂𝐣 to the interpolated electric field for
// the particle transport
   if(resistive_part) res = eta;

// Use the shape function to find the interpolated field. We think of particles
// "collecting" field from cells according to ShapeFunc().
   ShapeFunc(x / dx, Imax, i, f);

// magnetic field
   for(xyz = 1; xyz <= 3; xyz++) {
      Bint[xyz] = B[i[0]][xyz] * f[0];
      if(i[1]) Bint[xyz] += B[i[1]][xyz] * f[1];
      if(i[2]) Bint[xyz] += B[i[2]][xyz] * f[2];
   };

// electric field
   for(xyz = 1; xyz <= 3; xyz++) {
      Eint[xyz] = (E[i[0]][xyz] - res * j[i[0]][xyz])* f[0];
      if(i[1]) Eint[xyz] += (E[i[1]][xyz] - res * j[i[1]][xyz]) * f[1];
      if(i[2]) Eint[xyz] += (E[i[2]][xyz] - res * j[i[2]][xyz]) * f[2];
   };
};


//------------------------------------------------------------------------------
// Advance magnetic field using an implicit method
//------------------------------------------------------------------------------
// Input  | moments  | particle moments
// Input  | dt       | time step
//------------------------------------------------------------------------------

void fields_t::AdvanceImplicit(moments_t &moments, double dt)
{
   int i, xyz;
   double **a, **b, **c, **d, **x, mom[4], twodx, fourdx, dx2, Bndx2, dn4n;

   if(!isfld) return;

// Initialize matrix components
   a = new double *[Imax - 2];
   b = new double *[Imax - 2];
   c = new double *[Imax - 2];
   d = new double *[Imax - 2];
   x = new double *[Imax - 2];
   for(i = 0; i < Imax - 2; i++) {
      a[i] = new double [4];
      b[i] = new double [4];
      c[i] = new double [4];
      d[i] = new double [2];
      x[i] = new double [2];
   };

   twodx = 2.0 * dx;
   fourdx = 4.0 * dx;
   dx2 = dx * dx;

// Compute number density, bulk velocity, and electron pressure
   for(i = 2; i <= Imax - 1; i++) {
      moments.GetMoments(i, n[i], mom);
      for(xyz = 1; xyz <= 3; xyz++) u[i][xyz] = mom[xyz] / n[i];

// Adiabatic equation of state
      p[i] = pe0 * pow(n[i], gammaa);
   };

// Apply periodic conditions to moments
   n[1] = n[Imax - 1];
   for(xyz = 1; xyz <= 3; xyz++) u[1][xyz] = u[Imax - 1][xyz];
   p[1] = p[Imax - 1];
   n[Imax] = n[2];
   for(xyz = 1; xyz <= 3; xyz++) u[Imax][xyz] = u[2][xyz];
   p[Imax] = p[2];

// Compute the matrix
   for(i = 2; i <= Imax - 1; i++) {

      Bndx2 = B[i][1] / (dx2 * n[i]);
      dn4n = (n[i + 1] - n[i - 1]) / (4.0 * n[i]);

      a[i - 2][0] = -u[i - 1][1] / fourdx - eta / (2.0 * dx2);
      a[i - 2][1] = -Bndx2 / 2.0 * (1.0 + dn4n);
      a[i - 2][2] = -a[i - 2][1];
      a[i - 2][3] =  a[i - 2][0];
     
      b[i - 2][0] = 1.0 / dt + eta / dx2;
      b[i - 2][1] = Bndx2;
      b[i - 2][2] = -b[i - 2][1];
      b[i - 2][3] =  b[i - 2][0];

      c[i - 2][0] =  u[i + 1][1] / fourdx - eta / (2.0 * dx2);
      c[i - 2][1] = -Bndx2 / 2.0 * (1.0 - dn4n);
      c[i - 2][2] = -c[i - 2][1];
      c[i - 2][3] =  c[i - 2][0];

      d[i - 2][0] =  -a[i - 2][0] * B[i - 1][2] - a[i - 2][1] * B[i - 1][3]
                    - b[i - 2][0] * B[i][2]     - b[i - 2][1] * B[i][3]
                    - c[i - 2][0] * B[i + 1][2] - c[i - 2][1] * B[i + 1][3]
         + 2.0 / dt * B[i][2] + B[i][1] / twodx * (u[i + 1][2] - u[i - 1][2]);

      d[i - 2][1] =  -a[i - 2][2] * B[i - 1][2] - a[i - 2][3] * B[i - 1][3]
                    - b[i - 2][2] * B[i][2]     - b[i - 2][3] * B[i][3]          
                    - c[i - 2][2] * B[i + 1][2] - c[i - 2][3] * B[i + 1][3]
         + 2.0 / dt * B[i][3] + B[i][1] / twodx * (u[i + 1][3] - u[i - 1][3]);
   };

// Solve the matrix and assign values to B
   ThomasPeriodicBlock2(Imax - 2, a, b, c, d, x);
   for(i = 2; i <= Imax - 1; i++) {
      B[i][2] = x[i - 2][0];
      B[i][3] = x[i - 2][1];
   };

// Periodic boundaries for B
   B[1][2] = B[Imax - 1][2];
   B[1][3] = B[Imax - 1][3];
   B[Imax][2] = B[2][2];
   B[Imax][3] = B[2][3];

// Compute the current (this is used in "AdvanceCurrent()" to calculate E).
   for(i = 2; i <= Imax - 1; i++) {
      j[i][1] = 0.0;
      j[i][2] = (B[i - 1][3] - B[i + 1][3]) / twodx;
      j[i][3] = (B[i + 1][2] - B[i - 1][2]) / twodx;
   };

// Release memory
   for(i = 0; i < Imax - 2; i++) {
      delete[] a[i];
      delete[] b[i];
      delete[] c[i];
      delete[] d[i];
      delete[] x[i];
   };
   delete[] a;
   delete[] b;
   delete[] c;
   delete[] d;
   delete[] x;
};


//------------------------------------------------------------------------------
// Compute the electric field using Ohm's law and momentum advance
//------------------------------------------------------------------------------
// Input  | moments  | particle moments
// Input  | dt       | time step (typically 1/2 of the total)
//------------------------------------------------------------------------------

void fields_t::AdvanceElectric(moments_t &moments, double dt)
{
   int i, xyz;
   double mom[4], uB[4], jB[4], twodx;

   if(!isfld) return;
   twodx = 2.0 * dx;

// Re-compute number density, bulk velocity, and electron pressure for the
// free-streamings propagation
   for(i = 2; i <= Imax - 1; i++) {
      moments.GetMoments(i, n[i], mom);
      for(xyz = 1; xyz <= 3; xyz++) u[i][xyz] = mom[xyz] / n[i];

// Adiabatic equation of state
      p[i] = pe0 * pow(n[i], gammaa);
   };

// Apply periodic conditions to moments
   n[1] = n[Imax - 1];
   for(xyz = 1; xyz <= 3; xyz++) u[1][xyz] = u[Imax - 1][xyz];
   p[1] = p[Imax - 1];
   n[Imax] = n[2];
   for(xyz = 1; xyz <= 3; xyz++) u[Imax][xyz] = u[2][xyz];
   p[Imax] = p[2];

// Advance the ion momentum. Density and pressure are not changed because they
// are provided at the final position.
   for(i = 2; i <= Imax - 1; i++) {
      VectorProduct(j[i], B[i], jB);
      for(xyz = 1; xyz <= 3; xyz++) {
         u[i][xyz] += dt * (jB[xyz] - (p[i + 1] - p[i - 1]) / twodx) / n[i];
      };

      VectorProduct(u[i], B[i], uB);
      for(xyz = 1; xyz <= 3; xyz++) {
         E[i][xyz] = -uB[xyz] + jB[xyz] / n[i] + eta * j[i][xyz];
      };
      E[i][1] -= (p[i + 1] - p[i - 1]) / (twodx * n[i]);
   };

// Periodic boundaries for E
   for(xyz = 1; xyz <= 3; xyz++) {
      E[1][xyz] = E[Imax - 1][xyz];
      E[Imax][xyz] = E[2][xyz];
   };
};


//------------------------------------------------------------------------------
// Save the field arrays to a binary file
//------------------------------------------------------------------------------
// Input  | fname    | name of the file
// Input  | mark     | time mark for this output
//------------------------------------------------------------------------------

void fields_t::Save(const string &fname, double mark)
{
   ofstream datafile;

// Open new file for writing
   if(!isfld) return;
   datafile.open(fname.c_str(), ofstream::out | ofstream::binary);

// Write the time mark and the number of cells
   datafile.write((char *)&Imax, sizeof(int));
   datafile.write((char *)&mark, sizeof(double));

// Write out the field arrays
   datafile.write((char *)B[0], (Imax * 3 + 1) * sizeof(double));
   datafile.write((char *)E[0], (Imax * 3 + 1) * sizeof(double));

   datafile.close();
   cerr << "# Wrote field data\n";
};


//------------------------------------------------------------------------------
// Read the field arrays from a binary file
//------------------------------------------------------------------------------
// Input  | fname    | name of the file
// Output | mark     | time mark from the start of the file
// Return |          | 0 on success of 1 on failure
//------------------------------------------------------------------------------

int fields_t::Restore(const string &fname, double &mark)
{
   int imax_saved;
   ifstream datafile;

// Check the condition of the file
   if(!isfld) return 1;
   datafile.open(fname.c_str(), ifstream::in | ifstream::binary);
   if(!datafile.good()) {
      cerr << "# Error opening " << fname << endl;
      datafile.close();
      return 1;
   };

// Check that the number of cells in the file is the same
   datafile.read((char *)&imax_saved, sizeof(int));
   datafile.read((char *)&mark, sizeof(double));
   if(imax_saved != Imax) {
      cerr << "# Wrong cell count in " << fname << endl;
      datafile.close();
      return 1;
   };

// Read in the field arrays
   datafile.read((char *)B[0], (Imax * 3 + 1) * sizeof(double));
   datafile.read((char *)E[0], (Imax * 3 + 1) * sizeof(double));

   datafile.close();
   cerr << "# Read field data\n";
   
   return 0;
};


//------------------------------------------------------------------------------
// Print the fields as an ASCII file
//------------------------------------------------------------------------------
// Input  | scale    | scale factor for the fields
// Input  | fname    | name of the file
//------------------------------------------------------------------------------

void fields_t::Print(double scale, const string &fname)
{
   int i;

   ofstream textfile;
   if(!isfld) return;

// Open new file for writing
   textfile.open(fname.c_str(), ofstream::out);

   for(i = 2; i <= Imax - 1; i++) {
      textfile << setw(14) << setprecision(6) << (i - 1.5) * dx
               << setw(12) << setprecision(5) << B[i][1] * scale
               << setw(14) << setprecision(6) << B[i][2] * scale
               << setw(14) << setprecision(6) << B[i][3] * scale
               << setw(14) << setprecision(6) << E[i][1] * scale
               << setw(14) << setprecision(6) << E[i][2] * scale
               << setw(14) << setprecision(6) << E[i][3] * scale
               << endl;
   };

   textfile.close();
   cerr << "# Printed field data\n";
};


//------------------------------------------------------------------------------
// Print one field component for Fourier analysis
//------------------------------------------------------------------------------
// Input  | datafile | output file stream
// Input  | var      | tells which field component to store
//------------------------------------------------------------------------------

void fields_t::Dump(ofstream &datafile, int var)
{
   int i, xyz;
   double *field, **emf;
   
   if(!isfld || !datafile.good() || var < 1 || var > 6) return;
   field = new double[Imax - 2];

// find which component we need
   emf = (var <= 3 ? B : E);
   xyz = (var - 1) % 3 + 1;
   
// get one component from the EMF field array
   for(i = 2; i <= Imax - 1; i++) {
      field[i - 2] = emf[i][xyz];
   };
   datafile.write((char *)field, (Imax - 2) * sizeof(double));

   delete[] field;
};


//------------------------------------------------------------------------------
// Compute the energy in the EM fields
//------------------------------------------------------------------------------
// Return |          | total energy
//------------------------------------------------------------------------------

double fields_t::Energy(void)
{
   int i;
   double en = 0.0;

// Total energy is (𝐄²+𝐁²)/2
   for(i = 2; i <= Imax - 1; i++) {
      en += Norm2(E[i]) + Norm2(B[i]);
   };
   return en / 2.0;
};


//------------------------------------------------------------------------------
// Compute the fluctuating magnetic field variance
//------------------------------------------------------------------------------
// Return |          | field variance
//------------------------------------------------------------------------------

double fields_t::MagneticVariance(void)
{
   int i, xyz;
   double Bavg[4], delB2[4];
   
   memset(Bavg, 0, 4 * sizeof(double));
   memset(delB2, 0, 4 * sizeof(double));
   
// Compute the mean field
   for(i = 2; i <= Imax - 1; i++) {
      for(xyz = 1; xyz <= 3; xyz++) Bavg[xyz] += B[i][xyz] / (Imax - 2);
   };

// Compute the diagonal covariance tensor components
   for(i = 2; i <= Imax - 1; i++) {
      for(xyz = 1; xyz <= 3; xyz++) delB2[xyz] += Sqr(B[i][xyz] - Bavg[xyz]) / (Imax - 2);
   };

// The trace <dB_i*dB_i> is invariant under rotation
   return (delB2[1] + delB2[2] + delB2[3]) / Norm2(Bavg);
};


//------------------------------------------------------------------------------
// Compute the energy of the adiabatic electron population
//------------------------------------------------------------------------------
// Return |          | electron energy
//------------------------------------------------------------------------------

double fields_t::ElectronEnergy(void)
{
   int i;
   double en = 0.0;

// Electron energy is 𝒑/(𝛾-1)
   for(i = 2; i <= Imax - 1; i++) {
      en += p[i];
   };
   return en / (gammaa - 1.0);
};


//------------------------------------------------------------------------------
// Broadcast the fields to all nodes from the communicator
//------------------------------------------------------------------------------
// Input  | comm     | the communicator
// Input  | origin   | the origin (master or boss)
//------------------------------------------------------------------------------

void fields_t::Broadcast(int comm, int origin)
{
   MPI_Bcast(B[0], (Imax * 3 + 1), MPI_DOUBLE, origin, comm);
   MPI_Bcast(E[0], (Imax * 3 + 1), MPI_DOUBLE, origin, comm);
};


//------------------------------------------------------------------------------
// Class destructor - release field array memory
//------------------------------------------------------------------------------

fields_t::~fields_t()
{
   if(isfld) {
      Delete2D(B);
      Delete2D(E);
      Delete2D(j);
      delete[] n;
      Delete2D(u);
      delete[] p;
   };
};


//==============================================================================
// End of the file "hybrid_fields.cc"
//==============================================================================
