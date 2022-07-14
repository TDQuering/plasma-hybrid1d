//==============================================================================
// Start of the file "hybrid_distr.cc"
//
// Version 3d.
//
// Written by Vladimir Florinski
//==============================================================================


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include "geo_memory.hh"
#include "hybrid_distr.hh"

using namespace std;


//==============================================================================
// Standard analytic PDF expressions
//==============================================================================


// Uniform distribution, Cartesian
inline double PDF_range_0(double param1, double param2, double val)
{
   return (val > param1 && val <= param2 ? 1.0 : 0.0);
};


// Uniform distribution, cylindrical
inline double PDF_range_1(double param1, double param2, double val)
{
   return (val > param1 && val <= param2 ? val : 0.0);
};


// Uniform distribution, spherical
inline double PDF_range_2(double param1, double param2, double val)
{
   return (val > param1 && val <= param2 ? Sqr(val) : 0.0);
};


// Gaussian distribution, Cartesian
inline double PDF_Gauss_0(double param1, double param2, double val)
{
   return exp(-Sqr((val - param1) / param2));
};


// Gaussian distribution, cylindrical
inline double PDF_Gauss_1(double param1, double param2, double val)
{
   return val * exp(-Sqr((val - param1) / param2));
};


// Gaussian distribution, spherical
inline double PDF_Gauss_2(double param1, double param2, double val)
{
   return Sqr(val) * exp(-Sqr((val - param1) / param2));
};


//==============================================================================
// Custom distributions
//==============================================================================


const double in_ratio = 0.5; // ion to neutral ratio

// Pickup ion distribution without scattering
inline double PDF_pui_noscat(double param1, double param2, double val)
{
   double st, tanpsi;
   st = sqrt(1.0 - Sqr(val));
   tanpsi = param1 / sqrt(1.0 - Sqr(param1));
   return exp(param2 * ((1.0 - in_ratio) / st - tanpsi / val)) / val / st;
};


//==============================================================================
// "Dice roll" functions
//==============================================================================


// Constant distribution
inline double Roll_singular(double param1, double param2)
{
   return param1;
};


// Uniform distribution, Cartesian
inline double Roll_range_0(double param1, double param2)
{
   return param1 + drand48() * (param2 - param1);
};


// Uniform distribution, cylindrical
inline double Roll_range_1(double param1, double param2)
{
   return sqrt(Sqr(param1) + drand48() * (Sqr(param2) - Sqr(param1)));
};


// Uniform distribution, spherical
inline double Roll_range_2(double param1, double param2)
{
   return pow(Cube(param1) + drand48() * (Cube(param2) - Cube(param1)), 1.0 / 3.0);
};


// Gaussian distribution, Cartesian
inline double Roll_Gauss_0(double param1, double param2)
{
   return param1 + param2 * sqrt(-log(drand48())) * cos(2.0 * M_PI * drand48());
};


// Gaussian distribution, cylindrical
inline double Roll_Gauss_1(double param1, double param2)
{
   return param1 + param2 * sqrt(-log(drand48()));
};


// Gaussian distribution truncated to [-1,1]
inline double Roll_Gauss_R(double param1, double param2)
{
   double gauss_trial = -2.0;
   while(gauss_trial < -1.0 || gauss_trial > 1.0) {
      gauss_trial = param1 + param2 * sqrt(-log(drand48())) * cos(2.0 * M_PI * drand48());
   };
   return gauss_trial;
};


//==============================================================================
// The distribution_t class private methods
//==============================================================================


//------------------------------------------------------------------------------
// Binary search of the forward CDF table
//------------------------------------------------------------------------------
// Input  | table    | which table to use - first or second
// Input  | val      | random number giving the value of the CDF
//------------------------------------------------------------------------------

double distribution_t::SearchForwardCDF(int table, double val)
{
   int i1, i2, i3;
   double var_min, var_max, dvar, *cdf_table;
   
// determine which table we need
   if(table == 1) {
      var_min = limits[0];
      var_max = limits[1];
      cdf_table = cdf_table_f1;
   }
   else if(table == 2) {
      var_min = limits[2];
      var_max = limits[3];
      cdf_table = cdf_table_f2;
   }
   else {
      cerr << "# distribution_t: Invalid table request\n";
      return 0.0;
   };

// Bisection search (slow)
   dvar = (var_max - var_min) / table_size;
   i1 = 0;
   i2 = table_size;
   while(i2 - i1 > 1) {
      i3 = (i2 + i1) / 2;
      if(val < cdf_table[i3]) i2 = i3;
      else i1 = i3;
   };
   
   return var_min + (i1 * dvar * (cdf_table[i2] - val)
                  +  i2 * dvar * (val - cdf_table[i1]))
                  / (cdf_table[i2] - cdf_table[i1]);
};


//------------------------------------------------------------------------------
// Lookup from the inverse CDF table
//------------------------------------------------------------------------------
// Input  | table    | which table to use - first or second
// Input  | val      | random number giving the value of the CDF
//------------------------------------------------------------------------------

double distribution_t::SearchInverseCDF(int table, double val)
{
   int i;
   double dcdf, *cdf_table;
   
// determine which table we need
   if(table == 1) cdf_table = cdf_table_i1;
   else if(table == 2) cdf_table = cdf_table_i2;
   else {
      cerr << "# distribution_t: Invalid table request\n";
      return 0.0;
   };

// Linear interpolation (fast)
   dcdf = 1.0 / table_size;
   i = val / dcdf;
   return cdf_table[i] * (i + 1.0 - val / dcdf)
        + cdf_table[i + 1] * (val / dcdf - i);
};


//------------------------------------------------------------------------------
// Generates the forward CDF lookup tables
//------------------------------------------------------------------------------
// Input  | table    | which table to generate - first or second
//------------------------------------------------------------------------------

void distribution_t::GenerateForwardTables(int table)
{
   int i;
   double var, dvar, var_min, var_max, param1, param2, pdf_l, pdf_r;
   double *cdf_table, *norm;
   double (*PDF_func)(double param1, double param2, double val);

// determine which table we need
   if(table == 1) {
      var_min = limits[0];
      var_max = limits[1];
      param1 = params[0];
      param2 = params[1];
      norm = &norm_f1;
      PDF_func = PDF_func1;
   }
   else if(table == 2) {
      var_min = limits[2];
      var_max = limits[3];
      param1 = params[2];
      param2 = params[3];
      norm = &norm_f2;
      PDF_func = PDF_func2;
   }
   else {
      cerr << "# distribution_t: Invalid table request\n";
      return;
   };

// Singular PDFs can not be generated
   if(!PDF_func) return;

// Allocate memory for tables
   if(table == 1) {
      if(!cdf_table_f1) cdf_table_f1 = new double[table_size + 1];
      cdf_table = cdf_table_f1;
   }
   else {
      if(!cdf_table_f2) cdf_table_f2 = new double[table_size + 1];
      cdf_table = cdf_table_f2;
   };      

   dvar = (var_max - var_min) / table_size;
   cdf_table[0] = 0.0;
   pdf_r = PDF_func(param1, param2, var_min);

// Integrate the PDF
   for(i = 1; i <= table_size; i++) {
      var = var_min + i * dvar;
      pdf_l = pdf_r;
      pdf_r = PDF_func(param1, param2, var);
      cdf_table[i] = cdf_table[i - 1] + 0.5 * dvar * (pdf_l + pdf_r);
   };
   
// Normalize to unity
   *norm = cdf_table[table_size];
   for(i = 1; i <= table_size; i++) cdf_table[i] /= *norm;
};


//------------------------------------------------------------------------------
// Generates the inverse CDF lookup tables
//------------------------------------------------------------------------------
// Input  | table    | which table to generate - first or second
//------------------------------------------------------------------------------

void distribution_t::GenerateInverseTables(int table)
{
   int i, n1;
   double var, dvar, var_min, var_max, param1, param2, pdf_l, pdf_r, cdf, dcdf;
   double *cdf_table, *norm;
   double (*PDF_func)(double param1, double param2, double val);

// determine which table we need
   if(table == 1) {
      var_min = limits[0];
      var_max = limits[1];
      param1 = params[0];
      param2 = params[1];
      norm = &norm_i1;
      PDF_func = PDF_func1;
   }
   else if(table == 2) {
      var_min = limits[2];
      var_max = limits[3];
      param1 = params[2];
      param2 = params[3];
      norm = &norm_i2;
      PDF_func = PDF_func2;
   }
   else {
      cerr << "# distribution_t: Invalid table request\n";
      return;
   };

// Singular PDFs can not be generated
   if(!PDF_func) return;

// Allocate memory for tables
   if(table == 1) {
      if(!cdf_table_i1) cdf_table_i1 = new double[table_size + 1];
      cdf_table = cdf_table_i1;
   }
   else {
      if(!cdf_table_i2) cdf_table_i2 = new double[table_size + 1];
      cdf_table = cdf_table_i2;
   };      

// Compute the normalization
   dvar = (var_max - var_min) / table_size;
   pdf_r = PDF_func(param1, param2, var_min);
   *norm = 0.0;
   for(i = 1; i <= table_size; i++) {
      var = var_min + i * dvar;
      pdf_l = pdf_r;
      pdf_r = PDF_func(param1, param2, var);
      *norm += 0.5 * dvar * (pdf_l + pdf_r);
   };

// Integrate with a smaller step, recording values at regular intervals on the
// CDF axis
   n1 = table_size * icdf_steps;
   dcdf = 1.0 / table_size;
   pdf_r = PDF_func(param1, param2, var_min) / (*norm);
   dvar = (var_max - var_min) / n1;
   cdf = cdf_table[0] = 0.0;
   var = var_min;
   for(i = 1; i < table_size; i++) {
      while(cdf < i * dcdf) {
         var += dvar;
         pdf_l = pdf_r;
         pdf_r = PDF_func(param1, param2, var) / (*norm);
         cdf += 0.5 * dvar * (pdf_l + pdf_r);
      };
      cdf_table[i] = var;
   };
   cdf_table[table_size] = var_max;
};


//------------------------------------------------------------------------------
// Initialize a cylindrical distribution
//------------------------------------------------------------------------------
// Input  | p_inp[4] | parameters of the distribution
//------------------------------------------------------------------------------

void distribution_t::InitCyl(double *p_inp)
{
   memcpy(params, p_inp, 4 * sizeof(double));

// Singular v_perp - no PDF
   if(type1 == PDF_SINGULAR ) {
      limits[0] = limits[1] = params[0];
      PDF_func1 = NULL;
      Roll_func1 = Roll_singular;
   }

// Uniform v_perp - both PDF and roll function exist
   else if(type1 == PDF_RANGE) {
      limits[0] = params[0];
      limits[1] = params[1];
      PDF_func1 = PDF_range_1;
      Roll_func1 = Roll_range_1;
   }

// Gauss v_perp - no known roll function
   else if(type1 == PDF_GAUSS) {
      limits[0] = 0.0;
      limits[1] = params[0] + gauss_width * params[1];
      PDF_func1 = PDF_Gauss_1;
      Roll_func1 = NULL;
   }

// Centered Gaussian v_perp - both PDF and roll function exist
   else if(type1 == PDF_CGAUSS) {
      limits[0] = params[0] = 0.0;
      limits[1] = gauss_width * params[1];
      PDF_func1 = PDF_Gauss_1;
      Roll_func1 = Roll_Gauss_1;
   }

// Custom v_perp - not implemented
   else if(type1 == PDF_CUSTOM1) {
      limits[0] = 0.0;
      limits[1] = 0.0;
      PDF_func1 = NULL;
      Roll_func1 = NULL;
   }

// Either type zero or invalid input (beam default)
   else {
      type1 = PDF_ZERO;
      params[0] = params[1] = 0.0;
      limits[0] = limits[1] = 0.0;
      PDF_func1 = NULL;
      Roll_func1 = Roll_singular;
   };

//------------------------------------------------------------------------------
// Singular v_para - no PDF
   if(type2 == PDF_SINGULAR) {
      limits[2] = limits[3] = params[2];
      PDF_func2 = NULL;
      Roll_func2 = Roll_singular;
   }

// Uniform v_para - both PDF and roll function exist
   else if(type2 == PDF_RANGE) {
      limits[2] = params[2];
      limits[3] = params[3];
      PDF_func2 = PDF_range_0;
      Roll_func2 = Roll_range_0;
   }

// Gauss v_para - both PDF and roll function exist
   else if(type2 == PDF_GAUSS) {
      limits[2] = params[2] - gauss_width * params[3];
      limits[3] = params[2] + gauss_width * params[3];
      PDF_func2 = PDF_Gauss_0;
      Roll_func2 = Roll_Gauss_0;
   }

// Custom v_para - not implemented
   else if(type2 == PDF_CUSTOM1) {
      limits[2] = 0.0;
      limits[3] = 0.0;
      PDF_func2 = NULL;
      Roll_func2 = NULL;
   }

// Defaults to nondrifting annulus
   else {
      type2 = PDF_ZERO;
      params[2] = params[3] = 0.0;
      limits[2] = limits[3] = 0.0;
      PDF_func2 = NULL;
      Roll_func2 = Roll_singular;
   };
};


//------------------------------------------------------------------------------
// Initialize a spherical distribution
//------------------------------------------------------------------------------
// Input  | p_inp[4] | parameters of the distribution
//------------------------------------------------------------------------------

void distribution_t::InitSph(double *p_inp)
{
   memcpy(params, p_inp, 4 * sizeof(double));

// Singular v - no PDF
   if(type1 == PDF_SINGULAR) {
      limits[0] = limits[1] = params[0];
      PDF_func1 = NULL;
      Roll_func1 = Roll_singular;
   }

// Uniform v - both PDF and roll function exist
   else if(type1 == PDF_RANGE) {
      limits[0] = params[0];
      limits[1] = params[1];
      PDF_func1 = PDF_range_2;
      Roll_func1 = Roll_range_2;
   }

// Gauss v - roll function unavailable
   else if(type1 == PDF_GAUSS) {
      limits[0] = 0.0;
      limits[1] = params[0] + gauss_width * params[1];
      PDF_func1 = PDF_Gauss_2;
      Roll_func1 = NULL;
   }

// Defaults to stationary
   else {
      type1 = PDF_ZERO;
      params[0] = params[1] = 0.0;
      limits[0] = limits[1] = 0.0;
      PDF_func1 = NULL;
      Roll_func1 = Roll_singular;
   };

//------------------------------------------------------------------------------

// Singular mu - no PDF
   if(type2 == PDF_SINGULAR) {
      limits[2] = limits[3] = params[0];
      PDF_func2 = NULL;
      Roll_func2 = Roll_singular;
   }

// Uniform mu - both PDF and roll function exist
   else if(type2 == PDF_RANGE) {
      limits[2] = params[0];
      limits[3] = params[1];
      PDF_func2 = PDF_range_0;
      Roll_func2 = Roll_range_0;
   }

// Gauss mu - use truncated Gaussian roll
   else if(type2 == PDF_GAUSS) {
      limits[2] = -1.0;
      limits[3] = 1.0;
      PDF_func2 = PDF_Gauss_0;
      Roll_func2 = Roll_Gauss_R;
   }

// Custom mu (1)
   else if(type2 == PDF_CUSTOM1) {
      limits[2] = 1.0E-7;
      limits[3] = params[2];
      PDF_func2 = PDF_pui_noscat;
      Roll_func2 = NULL;
   }

// Defaults to nondrifting annulus
   else {
      type2 = PDF_ZERO;
      params[0] = params[1] = 0.0;
      limits[2] = limits[3] = 0.0;
      PDF_func2 = NULL;
      Roll_func2 = Roll_singular;
   };
};


//==============================================================================
// The distribution_t class public methods
//==============================================================================


//------------------------------------------------------------------------------
// Default constructor (cold stationary distribution)
//------------------------------------------------------------------------------

distribution_t::distribution_t()
{
   geom = false;
   type1 = 0;
   type2 = 1;
   memset(params, 0, n_par * sizeof(double));
   PDF_func1 = NULL;
   PDF_func2 = NULL;
   Roll_func1 = Roll_singular;
   Roll_func2 = Roll_singular;
   cdf_table_f1 = cdf_table_f2 = NULL;
   cdf_table_i1 = cdf_table_i2 = NULL;
};


//------------------------------------------------------------------------------
// Allocate memory, assign the PDF functions and compute the CDF tables
//------------------------------------------------------------------------------
// Input  | type     | type of distribution
// Input  | p_inp[4] | parameters of the distribution
//------------------------------------------------------------------------------

void distribution_t::Activate(int type, double *p_inp)
{
   geom = type / 100;
   type1 = (type - 100 * geom) / 10;
   type2 = (type - 100 * geom - 10 * type1);

// Initialize the parameter set
   if(geom == GEOM_CYL) InitCyl(p_inp);
   else if(geom == GEOM_SPH) InitSph(p_inp);
   else {
      cerr << "# distribution_t: Invalid geometry\n";
      return;
   };

// Generate the CDF tables   
   GenerateForwardTables(1);
   GenerateForwardTables(2);
   GenerateInverseTables(1);
   GenerateInverseTables(2);
};


//------------------------------------------------------------------------------
// Complete class constructor with initialization
//------------------------------------------------------------------------------
// Input  | type     | type of distribution
// Input  | p_inp[4] | parameters of the distribution
//------------------------------------------------------------------------------

distribution_t::distribution_t(int type, double *p_inp)
{
   Activate(type, p_inp);
};


//------------------------------------------------------------------------------
// Generates an instance of a three-dimensional velocity vector
//------------------------------------------------------------------------------
// Output | v[3]     | velocity vector
//------------------------------------------------------------------------------

void distribution_t::RollOne(double *v)
{
   double v1, v2, phi, rn1, rn2, st;

// stationary distribution as default
   memset(v, 0, 3 * sizeof(double));

// Select the best method to generate the distribution. Roll function is
// preferred, followed by inverse CDF, and finally forward CDF.
   rn1 = drand48();
   if(!Roll_func1 && !PDF_func1) return;
   else if(!Roll_func1 || (PDF_func1 && prefer_cdf)) {
      if(use_forward) v1 = SearchForwardCDF(1, rn1);
      else v1 = SearchInverseCDF(1, rn1);
   }
   else v1 = Roll_func1(params[0], params[1]);
   
   rn2 = drand48();
   if(!Roll_func2 && !PDF_func2) return;
   else if(!Roll_func2 || (PDF_func2 && prefer_cdf)) {
      if(use_forward) v2 = SearchForwardCDF(2, rn2);
      else v2 = SearchInverseCDF(2, rn2);
   }
   else v2 = Roll_func2(params[2], params[3]);

// Convert to Cartesian frame
   if(geom == GEOM_CYL) {
      if(type1 != PDF_ZERO) {
         phi = 2.0 * M_PI * drand48();
         v[0] = v1 * cos(phi);
         v[1] = v1 * sin(phi);
      };
      v[2] = v2;
   }
   else if(geom == GEOM_SPH && type1 != PDF_ZERO) {
      phi = 2.0 * M_PI * drand48();
      st = sqrt(1.0 - Sqr(v2));
      v[0] = v1 * st * cos(phi);
      v[1] = v1 * st * sin(phi);
      v[2] = v1 * v2;
   };
};


//------------------------------------------------------------------------------
// Class destructor
//------------------------------------------------------------------------------

distribution_t::~distribution_t()
{
   if(cdf_table_f1) delete[] cdf_table_f1;
   if(cdf_table_f2) delete[] cdf_table_f2;
   if(cdf_table_i1) delete[] cdf_table_i1;
   if(cdf_table_i2) delete[] cdf_table_i2;
};

