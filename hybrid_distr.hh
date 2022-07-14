//==============================================================================
// Start of the file "hybrid_distr.hh"
//
// Version 3d.
//
// Written by Vladimir Florinski
//==============================================================================


#ifndef HYBRID_DISTR
#define HYBRID_DISTR

#define GEOM_CYL     0
#define GEOM_SPH     1

#define PDF_ZERO     0
#define PDF_SINGULAR 1
#define PDF_RANGE    2
#define PDF_GAUSS    3
#define PDF_CGAUSS   4
#define PDF_CUSTOM1  5

#define is_cylindrical(x) (x / 100 == GEOM_CYL)
#define is_spherical(x)   (x / 100 == GEOM_SPH)


const int n_par = 4;                   // number of parameters per distribution
const int table_size = 10000;          // length of CDF lookup tables
const int icdf_steps = 100;            // CDF integration steps per dcdf
const double gauss_width = 3.0;        // width of a Gaussian (times sigma)
const bool prefer_cdf = false;         // true to prefer CDF lookup
const bool use_forward = true;         // true to use forward CDF lookup


//==============================================================================
// The distribution class
//==============================================================================


class distribution_t {

private:

   bool geom;             // geometry (cylindrical or spherical)
   int type1, type2;      // types (range, Gauss, etc.)
   double params[n_par];  // parameter list
   double limits[n_par];  // minimum and maximum values

// PDf and roll functions
   double (*PDF_func1) (double param1, double param2, double val);
   double (*PDF_func2) (double param1, double param2, double val);
   double (*Roll_func1)(double param1, double param2);
   double (*Roll_func2)(double param1, double param2);

// CDF tables
   double *cdf_table_f1;
   double *cdf_table_f2;
   double norm_f1, norm_f2;
   double *cdf_table_i1;
   double *cdf_table_i2;
   double norm_i1, norm_i2;

// Binary search of the forward CDF table
   double SearchForwardCDF(int table, double val);

// Lookup from the inverse CDF table
   double SearchInverseCDF(int table, double val);

// Initialize a cylindrical distribution
   void InitCyl(double *p_inp);

// Initialize a spherical distribution
   void InitSph(double *p_inp);

// Generates the forward CDF lookup tables
   void GenerateForwardTables(int table);

// Generates the inverse CDF lookup tables
   void GenerateInverseTables(int table);

public:

// Default constructor
   distribution_t();

// Allocate memory, assign the PDF functions and compute the CDF tables
   void Activate(int type, double *p_inp);

// Complete class constructor with initialization
   distribution_t(int type, double *p_inp);

// Generates an instance of a three-dimensional velocity vector
   void RollOne(double *v);

// Class destructor
   ~distribution_t();
};

#endif


//==============================================================================
// End of the file "hybrid_distr.hh"
//==============================================================================
