//==============================================================================
// Start of the file "geo_memory.hh"
//
// This file is part of the SPECTRUM suite version 5a.
//
// Written by Vladimir Florinski
//==============================================================================


#ifndef GEO_MEMORY
#define GEO_MEMORY
#include <cmath>

// macro definitions
#define MASTER 0
#define is_master !rank
#define LEFT 1
#define RIGHT 2

using namespace std;


// Common dimensionless constants
const double large = 1.0E20;
const double little = 1.0E-5;
const double small = 1.0E-8;
const double tiny = 1.0E-16;
const double pitwo = M_PI / 2.0;
const double twopi = 2.0 * M_PI;
const double fourpi = 4.0 * M_PI;
const double eightpi = 8.0 * M_PI;
const double sqrtpi = sqrt(M_PI);
const double sqrttwo = sqrt(2.0);

// 2D efficient array memory management
template <class T> T **Create2D(int n, int m)
{
   T **array = new T *[n + 1];
   array[0] = new T[n * m + 1]; array[1] = array[0];
   for(int i = 2; i <= n; i++) array[i] = array[i - 1] + m;
   return array;
};
template <class T> void Delete2D(T **array) {delete[] array[0]; delete[] array;};

template <class T> T **Map2D(int n, int m, T *start)
{
   T **array = new T *[n + 1];
   array[0] = start; array[1] = array[0];
   for(int i = 2; i <= n; i++) array[i] = array[i - 1] + m;
   return array;
};
template <class T> void Unmap2D(T **array) {delete[] array;};

// Inline routines for 2ⁿ, x², and x³
inline int Pow2(int n) {int p = 1; for(int i = 1; i <= n; i++, p *= 2); return p;};
inline double Sqr(double x) {return x * x;};
inline double Cube(double x) {return x * x * x;};

#endif


//==============================================================================
// End of the file "geo_memory.hh"
//==============================================================================
