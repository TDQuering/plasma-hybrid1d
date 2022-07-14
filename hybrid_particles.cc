//==============================================================================
// Start of the file "hybrid_particles.cc"
//
// Version 3b.
//
// Written by Vladimir Florinski
//==============================================================================


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <mpi.h>
#include "geo_memory.hh"
#include "geo_coord.hh"
#include "hybrid_particles.hh"

using namespace std;


//==============================================================================
// The particles_t class public methods
//==============================================================================


//------------------------------------------------------------------------------
// Allocate memory for particle coordinates
//------------------------------------------------------------------------------
// Input  | npart    | number of particles
// Input  | length   | size of the box
//------------------------------------------------------------------------------

void particles_t::Activate(int npart, double length)
{

// Can only be activated once
   if(isprt) return;

// Key parameters
   Lmax = npart;
   xmax = length;

// Memory for particle coordinates, velocities, and weights
   rvw = Create2D<double>(Lmax, 7);
   isprt = true;
};


//------------------------------------------------------------------------------
// Class constructor with memory allocation
//------------------------------------------------------------------------------
// Input  | npart    | number of particles
// Input  | length   | size of the box
//------------------------------------------------------------------------------
particles_t::particles_t(int npart, double length)
{
   isprt = ispdf = false;
   Activate(npart, length);
};


//------------------------------------------------------------------------------
// Load a gyrotropic particle distribution (cylindrical coordinates)
//------------------------------------------------------------------------------
// Input  | distr    | distribution class object
// Input  | d0       | number density, same as particle weight
// Input  | n0[4]    | symmetry axis
// Input  | cells    | number of cells (to generate uniform moments)
//------------------------------------------------------------------------------

void particles_t::Load(distribution_t &distr, double d0, double *n0, int cells)
{
   int l, ppc;
   double n1[4], n2[4], n3[4], vel[3];
   double dx, x_left;

   if(!isprt) return;

// Check if the number of particles is divisible by the number of cells
   ppc = Lmax / cells;
   if(Lmax % cells) {
      cerr << "# Unable to generate uniform moments with this Ncells\n";
      return;
   };
   dx = xmax / cells;

// Build a field-aligned velocity frame
   Copy(n0, n3);              // z' - along n0, in the XY plane
   Normalize(n3);
   VectorProduct(ez, n3, n1); // x' - perp. to n0, in the XY plane
   Normalize(n1);
   VectorProduct(n3, n1, n2); // y' - along z

// Generate uniform, gyrotropic PSD in the field-aligned system
   for(l = 1; l <= Lmax; l++) {

// Initially the number density is uniform. Integer division is used to obtain
// the left boundary of a cell.
      x_left = int((l - 1) / ppc) * dx;
      rvw[l][1] = x_left + dx * drand48();
      rvw[l][2] = 0.0;
      rvw[l][3] = 0.0;

// Request velocity roll from "distr"
      distr.RollOne(vel);

// Transform back to the fixed velocity frame
      rvw[l][4] = vel[0] * n1[1] + vel[1] * n2[1] + vel[2] * n3[1];
      rvw[l][5] = vel[0] * n1[2] + vel[1] * n2[2] + vel[2] * n3[2];
      rvw[l][6] = vel[0] * n1[3] + vel[1] * n2[3] + vel[2] * n3[3];

      rvw[l][7] = d0;
   };
};


//------------------------------------------------------------------------------
// Get a single particle
//------------------------------------------------------------------------------
// Input  | l        | particle's index
// Input  | r_part[4]| particle coordinate
// Input  | v_part[4]| particle velocity
// Input  | w_part   | particle weight
//------------------------------------------------------------------------------

void particles_t::GetOne(int l, double *r_part, double *v_part, double &w_part)
{
   if(l < 1 || l > Lmax || !isprt) {
      cerr << "# Invalid call to GetOne()\n";
      return;
   };

   Copy(&rvw[l][0], r_part);
   Copy(&rvw[l][3], v_part);
   w_part = rvw[l][7];
};


//------------------------------------------------------------------------------
// Generate a 3D Cartesian phase space density
//------------------------------------------------------------------------------
// Input  | size     | half-size of the PDF cube
// Input  | cells    | number of PDF cells in each dimension
// Input  | spc      | species number of the particle set
// file output written by Annaleah Ernst, summer 2015
//------------------------------------------------------------------------------

void particles_t::GeneratePDF(double size, int cells, int spc)
{
   int i, j, k, l, lrej;
   float dvel, incr;

   if(!isprt) return;

// If PDF was previously generated, release the memory because the new user
// may request a different PDF grid.
   if(ispdf) {
      for(i = 1; i <= Mmax; i++) Delete2D(PDF[i]);
      delete[] PDF;
      ispdf = false;
   };

// New PDF parameters
   vmax = size;
   Mmax = cells; // no need for ghost cells
   dvel = 2.0 * vmax / Mmax;
   incr = 1.0 / Lmax / Cube(dvel);

// We use the float type because precision is not an issue.
   PDF = new float **[Mmax + 1];
   for(i = 1; i <= Mmax; i++) {
      PDF[i] = Create2D<float>(Mmax, Mmax);
      memset(PDF[i][0], 0, (Mmax * Mmax + 1) * sizeof(float));
   };
   ispdf = true;
   lrej = 0;

// make a .bov file of the PDF
   ofstream hd_file, dat_file;
   string pdf_file_name = (spc < 10 ? "spc0" : "spc") + to_string(spc) + "pdf";

// generate header for a BOV file
   hd_file.open(pdf_file_name + ".bov", ofstream::out);
   if(!hd_file) {
      cerr << "Failed to open " << pdf_file_name << ".bov" << endl;
      return;
   }
   hd_file << "DATA_FILE: " << pdf_file_name << ".dat" << endl;
   hd_file << "DATA_SIZE: " << Mmax << " " << Mmax << " " << Mmax << endl;
   hd_file << "DATA_FORMAT: FLOAT" << endl;
   hd_file << "VARIABLE: particle" << endl;
   hd_file << "DATA_ENDIAN: LITTLE" << endl;
   hd_file << "CENTERING: zonal" << endl;
   hd_file << "BRICK_ORIGIN: " << -vmax << " " << -vmax << " " << -vmax << endl;
   hd_file << "BRICK_SIZE: " << 2.0 * vmax << " " << 2.0 * vmax << " " << 2.0 * vmax << endl;
   hd_file.close();

// Loop over particles
   for(l = 1; l <= Lmax; l++) {
      i = (rvw[l][4] + vmax) / dvel + 1;
      j = (rvw[l][5] + vmax) / dvel + 1;
      k = (rvw[l][6] + vmax) / dvel + 1;

// The distribution is normalized to 1
      if(i >= 1 && i <= Mmax && j >= 1 && j <= Mmax && k >= 1 && k <= Mmax) {
         PDF[i][j][k] += incr;
      }
      else lrej++;
   };

// open file
   dat_file.open(pdf_file_name + ".dat", ofstream::out | ofstream::binary);
   if(!dat_file) {
      cerr << "Failed to open " << pdf_file_name << ".dat" << endl;
      return;
   };

// write PDF to a BOV styled file
   for(i = 1; i <= Mmax; i++) {
      for(j = 1; j <= Mmax; j++) {
         for(k = 1; k <= Mmax; k++) {
            dat_file.write((char *)&PDF[i][j][k], sizeof(float));
         };
      };
   };
   dat_file.close();

// Print the number of particles that were outside the PDF cube
   cerr << "# Species " << spc << " is " << (double)lrej / (double)Lmax * 100.0
        << " %% outside the PDF cube\n";
};


//------------------------------------------------------------------------------
// Collect partial moments
//------------------------------------------------------------------------------
// InOut  | moments  | moments to modify
// Input  | update   | tells whether to rewrite the moments or add to existing
//------------------------------------------------------------------------------

void particles_t::ComputeMoments(moments_t &moments, bool update)
{
   int l;

   if(!isprt) return;

// Reset the moments if we are rewriting
   if(!update) moments.Reset();

// Add the particles to the moments one by one
   for(l = 1; l <= Lmax; l++) {
      moments.AddParticle(rvw[l][1], &rvw[l][3], rvw[l][7], true);
   };
};


//------------------------------------------------------------------------------
// Push particle velocity (leapfrog, matrix)
//------------------------------------------------------------------------------
// Input  | fields   | EM field
// Input  | dt       | time step (could be negative)
//------------------------------------------------------------------------------

void particles_t::PushMatrix(fields_t &fields, double dt)
{
   int l, cmp, cmq;
   double dt2, dt22, E[4], B[4], vB[4], D, M[4][4], R[4];

   if(!isprt) return;
   dt2 = dt / 2.0;
   dt22 = Sqr(dt2);

   for(l = 1; l <= Lmax; l++) {

// Interpolate fields to the particle's position
      fields.GetEMF(rvw[l][1], B, E);

// Integrate the Newton-Lorents equation. The velocity must be time-centered, so
// it is advanced from t-1/2 to t+1/2, whereas position is advanced from t to
// t+1. The force term vxB couples the velocity components. The matrix "M[][]"
// is the inverse of the N-L matrix.
      VectorProduct(&rvw[l][3], B, vB);
      M[1][1] = 1.0 + dt22 * Sqr(B[1]);
      M[1][2] = dt22 * B[1] * B[2] + dt2 * B[3];
      M[1][3] = dt22 * B[1] * B[3] - dt2 * B[2];
      M[2][1] = M[1][2] - dt * B[3];
      M[2][2] = 1.0 + dt22 * Sqr(B[2]);
      M[2][3] = dt22 * B[2] * B[3] + dt2 * B[1];
      M[3][1] = M[1][3] + dt * B[2];
      M[3][2] = M[2][3] - dt * B[1];
      M[3][3] = 1.0 + dt22 * Sqr(B[3]);

// The determinant
      D = M[1][1] + M[2][2] + M[3][3] - 2.0;

// The RHS (values at t-1/2)
      for(cmp = 1; cmp <= 3; cmp++) {
         R[cmp] = rvw[l][cmp + 3] + dt * E[cmp] + dt2 * vB[cmp];
      };

// Advance the velocity t-1/2 -> t+1/2: v=MR
      for(cmp = 1; cmp <= 3; cmp++) {
         rvw[l][cmp + 3] = 0.0;
         for(cmq = 1; cmq <= 3; cmq++) {
            rvw[l][cmp + 3] += M[cmp][cmq] * R[cmq] / D;
         };
      };
   };
};


//------------------------------------------------------------------------------
// Push particle velocity (leapfrog, Boris rotation)
//------------------------------------------------------------------------------
// Input  | fields   | EM field
// Input  | dt       | time step (could be negative)
//------------------------------------------------------------------------------

void particles_t::PushBoris(fields_t &fields, double dt)
{
   int l, cmp;
   double dt2, denom, E[4], B[4], vmin[4], vpls[4], vmt[4], vps[4], vprm[4],
          t[4], s[4];

   if(!isprt) return;
   dt2 = dt / 2.0;

   for(l = 1; l <= Lmax; l++) {

// Interpolate fields to the particle's position
      fields.GetEMF(rvw[l][1], B, E);
      
// Compute v- (half E push)
      for(cmp = 1; cmp <= 3; cmp++) {
         vmin[cmp] = rvw[l][cmp + 3] + E[cmp] * dt2;
      };

// Rotation
      for(cmp = 1; cmp <= 3; cmp++) t[cmp] = B[cmp] * dt2;
      denom = 2.0 / (1.0 + Norm2(t));
      for(cmp = 1; cmp <= 3; cmp++) s[cmp] = t[cmp] * denom;
      VectorProduct(vmin, t, vmt);
      for(cmp = 1; cmp <= 3; cmp++) vprm[cmp] = vmin[cmp] + vmt[cmp];
      VectorProduct(vprm, s, vps);
      for(cmp = 1; cmp <= 3; cmp++) vpls[cmp] = vmin[cmp] + vps[cmp];

// Compute v+ (half E push)
      for(cmp = 1; cmp <= 3; cmp++) {
         rvw[l][cmp + 3] = vpls[cmp] + E[cmp] * dt2;
      };
   };
};


//------------------------------------------------------------------------------
// Push particle velocity (Euler)
//------------------------------------------------------------------------------
// Input  | fields   | EM field
// Input  | dt       | time step (could be negative)
//------------------------------------------------------------------------------

void particles_t::PushV(fields_t &fields, double dt)
{
   int l, cmp;
   double E[4], B[4], vB[4];

   if(!isprt) return;

   for(l = 1; l <= Lmax; l++) {

// Interpolate fields to the particle's position
      fields.GetEMF(rvw[l][1], B, E);

// Advance the velocity via simple Euler scheme (should be used only to start
// a simulation).
      VectorProduct(&rvw[l][3], B, vB);
      for(cmp = 1; cmp <= 3; cmp++) {
         rvw[l][cmp + 3] += dt * (E[cmp] + vB[cmp]);
      };
   };
};


//------------------------------------------------------------------------------
// Push particle coordinate
//------------------------------------------------------------------------------
// Input  | dt       | time step (could be negative)
//------------------------------------------------------------------------------

void particles_t::PushX(double dt)
{
   int l, cmp;

   if(!isprt) return;

   for(l = 1; l <= Lmax; l++) {
      for(cmp = 1; cmp <= 3; cmp++) {
         rvw[l][cmp] += dt * rvw[l][cmp + 3];
      };

// periodic boundary in x
      if(rvw[l][1] <  0.0 ) rvw[l][1] += xmax;
      if(rvw[l][1] >= xmax) rvw[l][1] -= xmax;

// If the coordinate is still outside the box means particle traveled more than
// the box length in one step, which is clearly an error.
      if(rvw[l][1] < 0.0 || rvw[l][1] >= xmax) {
         cerr << "# Particle move error!\n";
      };
   };
};


//------------------------------------------------------------------------------
// Save the particle arrays to a binary file
//------------------------------------------------------------------------------
// Input  | datafile | output file stream
//------------------------------------------------------------------------------

void particles_t::Save(ofstream &datafile)
{
   if(!isprt || !datafile.good()) return;

   datafile.write((char *)&rvw[0][1], (Lmax * 7) * sizeof(double));
};


//------------------------------------------------------------------------------
// Read the particle arrays from a binary file
//------------------------------------------------------------------------------
// Input  | datafile | input file stream
//------------------------------------------------------------------------------

void particles_t::Restore(ifstream &datafile)
{
   if(!isprt || !datafile.good()) return;

// The binary data is a sequence of 7 numbers per particle. We need not be aware
// of CPU boundaries because the output is completely uniform.
   datafile.read((char *)&rvw[0][1], (Lmax * 7) * sizeof(double));
};


//------------------------------------------------------------------------------
// Print particle locations as projections onto a phase space plane
//------------------------------------------------------------------------------
// Input  | type     | type of desired output (1, 2, 3)
// Input  | skip     | skip this many particles per output
// Input  | scale    | scale factor for velocity
// Input  | n3[4]    | unit vector in the direction of velocity coordinate
//        |          | or normal to the plane of output
// Input  | fname    | name of the file
//------------------------------------------------------------------------------

void particles_t::Print(int type, int skip, double scale, const double *n3, const string &fname)
{
   int l, count = 0;
   double vel_1, vel_2, rad_1, rad_2, n1[4], n2[4];
   ofstream textfile;

   if(!isprt) return;

// Compute the normal plane
   if(fabs(n3[1]) > fmax(fabs(n3[2]), fabs(n3[3]))) VectorProduct(n3, ey, n2);
   else if(fabs(n3[2]) > fabs(n3[3])) VectorProduct(n3, ez, n2);
   else VectorProduct(n3, ex, n2);
   Normalize(n2);
   VectorProduct(n2, n3, n1);

// Open new file for writing
   textfile.open(fname.c_str(), ofstream::out);

   for(l = 1; l <= Lmax; l += skip) {
   
// r1-r2 plot - project coordinate onto the normal plane
      if(type == 1) {
         rad_1 = ScalarProduct(&rvw[l][0], n1);
         rad_2 = ScalarProduct(&rvw[l][0], n2);
         textfile << setw(14) << setprecision(6) << rad_1
                  << setw(14) << setprecision(6) << rad_2
                  << endl;
      }

// x-v plot - project velocity onto the normal
      if(type == 2) {
         vel_1 = ScalarProduct(&rvw[l][3], n3);
         textfile << setw(14) << setprecision(6) << rvw[l][1]
                  << setw(14) << setprecision(6) << vel_1 * scale
                  << endl;
      }

// v1-v2 plot - project velocity onto the normal plane
      else {
         vel_1 = ScalarProduct(&rvw[l][3], n1);
         vel_2 = ScalarProduct(&rvw[l][3], n2);
         textfile << setw(14) << setprecision(6) << vel_1 * scale
                  << setw(14) << setprecision(6) << vel_2 * scale
                  << endl;
      };
      count++;
   };

   textfile.close();
   cerr << "# Printed a map with " << count << " particles\n";
};


//------------------------------------------------------------------------------
// Compute the kinetic energy of all particles
//------------------------------------------------------------------------------
// Return |          | total energy
//------------------------------------------------------------------------------

double particles_t::Energy(void)
{
   int l;
   double en = 0.0;

   for(l = 1; l <= Lmax; l++) {
      en += Norm2(&rvw[l][3]) * rvw[l][7];
   };
   return en / 2.0;
};


//------------------------------------------------------------------------------
// Compute the mean and variance of the pitch angle
//------------------------------------------------------------------------------
// Input  | n0[4]    | z-axis unit vector
// Output | mu[4]    | average pitch angle cosine <mu>
// Output | mu2[4]   | average square of the pitch angle cosine <mu^2>
//------------------------------------------------------------------------------

void particles_t::MuAverage(double *n0, double &mu, double &mu2)
{
   int l;
   double vel_n;

   mu = mu2 = 0.0;
   for(l = 1; l <= Lmax; l++) {
      vel_n = Norm(&rvw[l][3]);
      mu += ScalarProduct(&rvw[l][3], n0) / vel_n * rvw[l][7];
      mu2 += Sqr(ScalarProduct(&rvw[l][3], n0) / vel_n) * rvw[l][7];
   };
};


//------------------------------------------------------------------------------
// Send the particle data to another CPU
//------------------------------------------------------------------------------
// Input  | dest     | destination rank
//------------------------------------------------------------------------------

void particles_t::Send(int dest)
{
   MPI_Send(rvw[0], Lmax * 7 + 1, MPI_DOUBLE, dest, 101, MPI_COMM_WORLD);
};


//------------------------------------------------------------------------------
// Receive the particle data from another CPU
//------------------------------------------------------------------------------
// Input  | source   | source rank
//------------------------------------------------------------------------------

void particles_t::Receive(int source)
{
   MPI_Status status;

   MPI_Recv(rvw[0], Lmax * 7 + 1, MPI_DOUBLE, source, 101, MPI_COMM_WORLD, &status);
};


//------------------------------------------------------------------------------
// Copy particle coordinates, velocities, and weights
//------------------------------------------------------------------------------
// Input  | origin   | particle structure to copy
// Return |          | this with a copy of origin
//------------------------------------------------------------------------------

particles_t &particles_t::operator =(particles_t &origin)
{
   if(isprt && origin.isprt && (Lmax == origin.Lmax)) {
      memcpy(rvw[0], origin.rvw[0], (Lmax * 7 + 1) * sizeof(double));
   };

   return *this;
};


//------------------------------------------------------------------------------
// Class destructor - release particle array memory
//------------------------------------------------------------------------------

particles_t::~particles_t()
{
   int i;

   if(isprt) Delete2D(rvw);

// PDF array
   if(ispdf) {
      for(i = 1; i <= Mmax; i++) Delete2D(PDF[i]);
      delete[] PDF;
   };
};


//==============================================================================
// End of the file "hybrid_particles.cc"
//==============================================================================
