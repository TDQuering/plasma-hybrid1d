//==============================================================================
// Start of the file "hybrid_fields.hh"
//
// Version 3d.
//
// Written by Vladimir Florinski
//==============================================================================


#ifndef HYBRID_FIELDS
#define HYBRID_FIELDS

#include "hybrid_moments.hh"
#include <complex>


// adiabatic index
const double gammaa = 5.0 / 3.0;
//const double gammaa = 1.001; // isothermal

// resistivity: 10^-6 is a reasonable number
const double eta = 0.0;

// whether to add resistive part to E for as used by the particles
const bool resistive_part = false;


//==============================================================================
// The fields class
//==============================================================================

class fields_t {

private:

   int    Imax;   // number of grid cells
   double xmax;   // size of the box
   double dx;     // size of a grid cell
   double pe0;    // initial electron pressure
   double wpiwci; // plasma to cyclotron frequency ratio
   bool   isfld;  // flag telling whether field arrays were generated

   double **B;    // magnetic field
   double **E;    // electric field
   double **j;    // current density

   double *n;     // number density
   double **u;    // bulk velocity
   double *p;     // electron pressure

public:

// Get the number of cells
   int GetCells(void) {return Imax;};

// Get the length of the box
   int GetLength(void) {return xmax;};

// Default constructor
   fields_t() {isfld = false;};

// Allocate memory for field arrays
   void Activate(int cells, double length, double fratio);

// Class constructor with memory allocation
   fields_t(int cells, double length, double fratio);

// Apply the initial condition for the EM field
   void Init(int type, double *B0, double epre);

// Return the interpolated EM fields
   void GetEMF(double x, double *Bint, double *Eint);

// Advance magnetic field using an implicit method
   void AdvanceImplicit(moments_t &moments, double dt);

// Apply the driving current to produce the magnetosonic wave
   void ApplyDrivingCurrent(double t, double omega_t, double j_amp);

// Advance the fields according to a wave, decoupled from the moments
   void fields_t::AdvanceDecoupled(double t, double wave_omega, double wave_amp, std::complex<double> Fluct_E[4], double B0[4], double mu0);

// Compute the electric field using Ohm's law and momentum advance
   void AdvanceElectric(moments_t &moments, double dt);

// Save the field arrays to a binary file
   void Save(const string &fname, double mark);

// Read the field arrays from a binary file
   int Restore(const string &fname, double &mark);

// Print the fields as an ASCII file
   void Print(double scale, const string &fname);

// Print one field component for Fourier analysis
   void Dump(ofstream &datafile, int var);

// Compute the energy in the EM fields
   double Energy(void);

// Compute the fluctuating magnetic field variance
   double MagneticVariance(void);

// Compute the energy of the adiabatic electron population
   double ElectronEnergy(void);

// Broadcast the fields to all members of the given communicator
   void Broadcast(int comm, int origin);

// Class destructor - release field array memory
   ~fields_t();
};

#endif


//==============================================================================
// End of the file "hybrid_fields.hh"
//==============================================================================
