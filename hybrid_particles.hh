//==============================================================================
// Start of the file "hybrid_particles.hh"
//
// Version 3d.
//
// Written by Vladimir Florinski
//==============================================================================


#ifndef HYBRID_PARTICLES
#define HYBRID_PARTICLES

#include "hybrid_distr.hh"
#include "hybrid_fields.hh"
#include "hybrid_moments.hh"


// unit vectors
const double ex[4] = {0.0, 1.0, 0.0, 0.0};
const double ey[4] = {0.0, 0.0, 1.0, 0.0};
const double ez[4] = {0.0, 0.0, 0.0, 1.0};


//==============================================================================
// The particles class
//==============================================================================

class particles_t {

private:

   int    Lmax;   // number of particles
   double xmax;   // size of the box (only meaningful for periodic boundaries)

   double **rvw;  // particle coordinate, velocity, and weight
   bool   isprt;  // flag telling whether particle arrays were generated
   
   int    Mmax;   // number of PDF cells in each dimension
   double vmax;   // half-size of the PDF cube
   float  ***PDF; // phase space density array
   bool   ispdf;  // flag telling whether PDF array was generated

public:

// Return the number of particles
   int GetNpart(void) {return Lmax;};

// Get the length of the box
   int GetLength(void) {return xmax;};

// Default constructor
   particles_t() {isprt = ispdf = false;};

// Allocate memory for particle coordinates
   void Activate(int npart, double length);

// Class constructor with memory allocation
   particles_t(int npart, double length);

// Load a gyrotropic velocity distribution
   void Load(distribution_t &distr, double d0, double *n0, int cells);

// Get a single particle
   void GetOne(int l, double *r_part, double *v_part, double &w_part);

// Generate a 3D Cartesian phase space density
   void GeneratePDF(double size, int cells, int spc);

// Collect partial moments
   void ComputeMoments(moments_t &moments, bool update);

// Push particle velocity (leapfrog, matrix)
   void PushMatrix(fields_t &fields, double dt);

// Push particle velocity (leapfrog, Boris rotation)
   void PushBoris(fields_t &fields, double dt);

// Push particle velocity (Euler)
   void PushV(fields_t &fields, double dt);

// Push particle coordinate
   void PushX(double dt);

// Save the particle arrays to a binary file
   void Save(ofstream &datafile);

// Read the particle arrays from a binary file
   void Restore(ifstream &datafile);

// Print particle locations as projections onto a phase space plane
   void Print(int type, int skip, double scale, const double *n3, const string &fname);

// Compute the kinetic energy of all particles
   double Energy(void);

// Compute the mean and variance of the pitch angle
   void MuAverage(double *n0, double &mu, double &mu2);

// Send the particle data to another CPU
   void Send(int dest);

// Receive the particle data from another CPU
   void Receive(int source);

// Copy particle coordinates, velocities, and weights
   particles_t &operator =(particles_t &origin);

// Class destructor - release particle array memory
   ~particles_t();
};

#endif


//==============================================================================
// End of the file "hybrid_particles.hh"
//==============================================================================
