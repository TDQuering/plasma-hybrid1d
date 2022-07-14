//==============================================================================
// Start of the file "hybrid_moments.hh"
//
// Version 3b.
//
// Written by Vladimir Florinski
//==============================================================================


#ifndef HYBRID_MOMENTS
#define HYBRID_MOMENTS

// shape function order
const int shpfc = 2;

// Shape function for particle-to-field and field-to-particle mappings
void ShapeFunc(double xrel, int n, int *i, double *f);


//==============================================================================
// The moments class
//==============================================================================

class moments_t {

private:

   int    Imax;   // number of grid cells
   double xmax;   // size of the box
   double dx;     // size of a grid cell
   bool   ismom;  // flag telling whether moment arrays were generated

// moments
   double  *den;  // number density
   double **mom;  // momentum
   double **enr;  // energy

public:

// Get the number of cells
   int GetCells(void) {return Imax;};

// Get the length of the box
   int GetLength(void) {return xmax;};

// Default constructor
   moments_t() {ismom = false;};
   
// Allocate memory for moments arrays
   void Activate(int cells, double length);

// Class constructor with memory allocation
   moments_t(int cells, double length);

// Reset all moments to zero
   void Reset(void);

// Add one particle to the moment arrays
   void AddParticle(double x, double *v, double wgt, bool doenrgy);

// Apply a 3-point digital fitler
   void Filter(double w);

// Return density and velocity in a cell
   void GetMoments(int i, double &denc, double *momc);

// Return total grid mass, momentum, and energy
   void GetTotal(double &mass, double *momt, double &enrg);

// Print the moments as an ASCII file
   void Print(double scale, const string &fname);

// Collect all moments
   void Collect(int comm, int root);

// Add partial moments to the current destribution
   moments_t &operator +=(moments_t &momentp);

// Average two sets of moments with given weight
   void Average(moments_t &moment1, moments_t &moment2, double weight1);

// Class destructor - release moments array memory
   ~moments_t();
};

#endif


//==============================================================================
// End of the file "hybrid_moments.hh"
//==============================================================================
