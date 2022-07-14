//==============================================================================
// Start of the file "hybrid_moments.cc"
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
#include "hybrid_moments.hh"

using namespace std;


//------------------------------------------------------------------------------
// Shape function for particle-to-field and field-to-particle mappings
//------------------------------------------------------------------------------
// Input  | xrel     | coordinate in units of dx
// Input  | n        | total number of cells (including ghost)
// Output | i[3]     | cells to be updated
// Output | f[3]     | fraction going into each cell
//------------------------------------------------------------------------------

void ShapeFunc(double xrel, int n, int *i, double *f)
{
// The point x=xmax (exactly) is mapped to x=0, so we use the ">=" operation
// for the right boundary. This simplifies the following code computing i and f.
   if(xrel < 0.0 || xrel >= n - 2) {
      cerr << "# Shape function error!\n";
   };

// The weight is distributed over 3 cells: i₀, i₁, i₂. For cell centers we have
// x(i₀)=i₀-3/2, x(i₁)=i₀-1/2, x(i₂)=i₀+1/2
   switch(shpfc) {

// cloud in cell (1 order) interpolation
   case 1:
      i[0] = xrel + 1.5;
      i[1] = i[0] + 1;
      i[2] = 0;
      f[0] = i[0] - xrel - 0.5;
      f[1] = 1.0 - f[0];
      f[2] = 0.0;

      if(i[0] == 1) i[0] = n - 1;
      if(i[1] == n) i[1] = 2;
      break;

// triangular shape cloud (2 order) interpolation
   case 2:
      i[0] = xrel + 1.0;
      i[1] = i[0] + 1;
      i[2] = i[1] + 1;
      f[0] = Sqr(i[0] - xrel) / 2.0;
      f[2] = Sqr(i[0] - xrel - 1.0) / 2.0;
      f[1] = 1.0 - f[0] - f[2];

      if(i[0] == 1) i[0] = n - 1;
      if(i[2] == n) i[2] = 2;
      break;
      
// nearest grid point (0 order) interpolation
   default:
      i[0] = xrel + 2.0;
      i[1] = 0;
      i[2] = 0;
      f[0] = 1.0;
      f[1] = 0.0;
      f[2] = 0.0;
      break;
   };
};


//==============================================================================
// The moments_t class public methods
//==============================================================================


//------------------------------------------------------------------------------
// Allocate memory for moments arrays
//------------------------------------------------------------------------------
// Input  | cells    | number of grid cells
// Input  | length   | size of the box
//------------------------------------------------------------------------------

void moments_t::Activate(int cells, double length)
{

// Can only be activated once
   if(ismom) return;

// Key parameters
   Imax = cells + 2; // 2 ghost cells
   xmax = length;
   dx = xmax / (Imax - 2);

// Memory for density, momentum, and energy arrays
   den = new double[Imax + 1];
   mom = Create2D<double>(Imax, 3);
   enr = Create2D<double>(Imax, 3);

// Initialize to zero
   Reset();
   ismom = true;
};


//------------------------------------------------------------------------------
// Class constructor with memory allocation
//------------------------------------------------------------------------------
// Input  | cells    | number of grid cells
// Input  | length   | size of the box
//------------------------------------------------------------------------------

moments_t::moments_t(int cells, double length)
{
   ismom = false;
   Activate(cells, length);
};


//------------------------------------------------------------------------------
// Reset all moments to zero
//------------------------------------------------------------------------------

void moments_t::Reset(void)
{
   if(!ismom) return;

   memset(den, 0, (Imax + 1) * sizeof(double));
   memset(mom[0], 0, (3 * Imax + 1) * sizeof(double));
   memset(enr[0], 0, (3 * Imax + 1) * sizeof(double));
};


//------------------------------------------------------------------------------
// Add one particle to the moment arrays
//------------------------------------------------------------------------------
// Input  | x        | particle coordinate
// Input  | v        | particle velocity
// Input  | wgt      | particle weight
// Input  | doenrgy  | request to compute the energy
//------------------------------------------------------------------------------

void moments_t::AddParticle(double x, double *v, double wgt, bool doenrgy)
{
   int cmp, i[3];
   double moment, f[3];
   if(!ismom) return;

// Distribute particle's weight according to the shape function (which
// guarantees correct treatment of the boundaries).
   ShapeFunc(x / dx, Imax, i, f);

// density
   den[i[0]] += wgt * f[0];
   if(i[1]) den[i[1]] += wgt * f[1];
   if(i[2]) den[i[2]] += wgt * f[2];

// momentum
   for(cmp = 1; cmp <= 3; cmp++) {
      moment = v[cmp] * wgt;
      mom[i[0]][cmp] += moment * f[0];
      if(i[1]) mom[i[1]][cmp] += moment * f[1];
      if(i[2]) mom[i[2]][cmp] += moment * f[2];
   };

// energy, only if requested
   for(cmp = 1; doenrgy && cmp <= 3; cmp++) {
      moment = Sqr(v[cmp]) * wgt;
      enr[i[0]][cmp] += moment * f[0];
      if(i[1]) enr[i[1]][cmp] += moment * f[1];
      if(i[2]) enr[i[2]][cmp] += moment * f[2];
   };
};


//------------------------------------------------------------------------------
// Apply a 3-point digital fitler
//------------------------------------------------------------------------------
// Input  | w         | weight
//------------------------------------------------------------------------------

void moments_t::Filter(double w)
{
   int i, cmp;
   double *den0, **men0, denom;
   if(!ismom) return;
   denom = 1.0 + 2.0 * w;

   den0 = new double[Imax + 1];
   men0 = Create2D<double>(Imax, 3);

// apply periodic boundary conditions
   den[1] = den[Imax - 1];
   den[Imax] = den[2];
   for(cmp = 1; cmp <= 3; cmp++) {
      mom[1][cmp] = mom[Imax - 1][cmp];
      mom[Imax][cmp] = mom[2][cmp];
      enr[1][cmp] = enr[Imax - 1][cmp];
      enr[Imax][cmp] = enr[2][cmp];
   };

// filter the density
   memcpy(den0, den, (Imax + 1) * sizeof(double));
   for(i = 2; i <= Imax - 1; i++) {
      den[i] = (den0[i] + w * (den0[i - 1] + den0[i + 1])) / denom;
   };

// filter the momentum
   memcpy(men0[0], mom[0], (3 * Imax + 1) * sizeof(double));
   for(i = 2; i <= Imax - 1; i++) {
      for(cmp = 1; cmp <= 3; cmp++) {
         mom[i][cmp] = (men0[i][cmp] + w * (men0[i - 1][cmp] + men0[i + 1][cmp])) / denom;
      };
   };
   
// energy filtering is optional (for diagnostic only)
   memcpy(men0[0], enr[0], (3 * Imax + 1) * sizeof(double));
   for(i = 2; i <= Imax - 1; i++) {
      for(cmp = 1; cmp <= 3; cmp++) {
         enr[i][cmp] = (men0[i][cmp] + w * (men0[i - 1][cmp] + men0[i + 1][cmp])) / denom;
      };
   };

   delete[] den0;
   Delete2D(men0);
};


//------------------------------------------------------------------------------
// Return density and velocity in a cell
//------------------------------------------------------------------------------
// Input  | i         | cell index
// Output | denc      | number density
// Output | momc[4]   | momentum
//------------------------------------------------------------------------------

void moments_t::GetMoments(int i, double &denc, double *momc)
{
   denc = den[i];
   memcpy(momc, mom[i], 4 * sizeof(double));
};


//------------------------------------------------------------------------------
// Return total grid mass, momentum, and energy
//------------------------------------------------------------------------------
// Output | mass      | mass
// Output | momt[4]   | momentum
// Output | enrg      | energy
//------------------------------------------------------------------------------

void moments_t::GetTotal(double &mass, double *momt, double &enrg)
{
   int i, cmp;
   mass = enrg = 0.0;
   memset(momt, 0, 4 * sizeof(double));
   
   for(i = 2; i <= Imax - 1; i++) {
      mass += den[i];
      for(cmp = 1; cmp <= 3; cmp++) {
         momt[cmp] += mom[i][cmp];
         enrg += enr[i][cmp];
      };
   };
};


//------------------------------------------------------------------------------
// Print the moments as an ASCII file
//------------------------------------------------------------------------------
// Input  | scale    | scale factor for velocity
// Input  | fname    | name of the file
//------------------------------------------------------------------------------

void moments_t::Print(double scale, const string &fname)
{
   int i, cmp;
   double v_avg[4], v2_avg[4];
   ofstream textfile;
   if(!ismom) return;

// Open new file for writing
   textfile.open(fname.c_str(), ofstream::out);

// Print the number density, bulk velocity, and random speed
   for(i = 2; i <= Imax - 1; i++) {
      for(cmp = 1; cmp <= 3; cmp++) {
         v_avg[cmp] = mom[i][cmp] / den[i];
         v2_avg[cmp] = enr[i][cmp] / den[i] - Sqr(v_avg[cmp]);
      };
      textfile << setw(14) << setprecision(6) << (i - 1.5) * dx
               << setw(14) << setprecision(6) << den[i]
               << setw(14) << setprecision(6) << v_avg[1] * scale
               << setw(14) << setprecision(6) << v_avg[2] * scale
               << setw(14) << setprecision(6) << v_avg[3] * scale
               << setw(14) << setprecision(6) << sqrt(v2_avg[1]) * scale
               << setw(14) << setprecision(6) << sqrt(v2_avg[2]) * scale
               << setw(14) << setprecision(6) << sqrt(v2_avg[3]) * scale
               << endl;
   };

   textfile.close();
};


//------------------------------------------------------------------------------
// Collect all moments
//------------------------------------------------------------------------------
// Input  | comm     | the communicator
// Input  | root     | the destination (master or boss)
//------------------------------------------------------------------------------
void moments_t::Collect(int comm, int root)
{
   int rank;

   MPI_Comm_rank(comm, &rank);
   MPI_Reduce((rank == root ? MPI_IN_PLACE : den), den, Imax + 1, MPI_DOUBLE, MPI_SUM, root, comm);
   MPI_Reduce((rank == root ? MPI_IN_PLACE : mom[0]), mom[0], (Imax * 3 + 1), MPI_DOUBLE, MPI_SUM, root, comm);
   MPI_Reduce((rank == root ? MPI_IN_PLACE : enr[0]), enr[0], (Imax * 3 + 1), MPI_DOUBLE, MPI_SUM, root, comm);
};


//------------------------------------------------------------------------------
// Add partial moments to the current destribution
//------------------------------------------------------------------------------
// Input  | momentp  | partial moments to add
// Return |          | this incremented with "momentp"
//------------------------------------------------------------------------------

moments_t &moments_t::operator +=(moments_t &momentp)
{
   int i, cmp;

   if(ismom && momentp.ismom) {
      for(i = 2; i <= Imax - 1; i++) {
         den[i] += momentp.den[i];
         for(cmp = 1; cmp <= 3; cmp++) {
            mom[i][cmp] += momentp.mom[i][cmp];
            enr[i][cmp] += momentp.enr[i][cmp];
         };
      };
   };

   return *this;
};


//------------------------------------------------------------------------------
// Average two sets of moments with given weight
//------------------------------------------------------------------------------
// Input  | moment1  | first set of moments
// Input  | moment2  | second set of moments
// Input  | weight1  | Weight for the first moment
//------------------------------------------------------------------------------

void moments_t::Average(moments_t &moment1, moments_t &moment2, double weight1)
{
   int i, cmp;
   double weight2;
   weight2 = 1.0 - weight1;
   if(!ismom || !moment1.ismom || !moment2.ismom) return;

   for(i = 2; i <= Imax - 1; i++) {
      den[i] = weight1 * moment1.den[i] + weight2 * moment2.den[i];
      for(cmp = 1; cmp <= 3; cmp++) {
         mom[i][cmp] = weight1 * moment1.mom[i][cmp] + weight2 * moment2.mom[i][cmp];
         enr[i][cmp] = weight1 * moment1.enr[i][cmp] + weight2 * moment2.enr[i][cmp];
      };
   };
};


//------------------------------------------------------------------------------
// Class destructor - release moments array memory
//------------------------------------------------------------------------------

moments_t::~moments_t()
{
   if(ismom) {
      delete[] den;
      Delete2D(mom);
      Delete2D(enr);
   };
};


//==============================================================================
// End of the file "hybrid_moments.cc"
//==============================================================================
