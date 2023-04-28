//==============================================================================
// Start of the file "geo_coord.hh"
//
// This file is part of the SPECTRUM suite version 5a.
//
// Written by Vladimir Florinski
//==============================================================================


#ifndef GEO_COORD
#define GEO_COORD

#include <cmath>
#include <cstring>
#include <complex>

using namespace std;


// Convention for spherical coordinates: latitude = 𝜋/2 - 𝜃, longitude = 𝜑

const double radial_unit_vec[4] = {0.0, 1.0, 0.0, 0.0};


// Convert from radians to degrees
inline double RadToDeg(double rad) {return rad * 180.0 / M_PI;};

// Convert from degrees to radians
inline double DegToRad(double deg) {return deg / 180.0 * M_PI;};


//==============================================================================
// Cartesian-polar coordinate transformation
//==============================================================================


// Convert a vector from polar to Cartesian (return x coordinate)
double GetUx(double ur, double ut, double up, double sintheta, double costheta,
   double sinphi, double cosphi);

// Convert a vector from polar to Cartesian (return y coordinate)
double GetUy(double ur, double ut, double up, double sintheta, double costheta,
   double sinphi, double cosphi);

// Convert a vector from polar to Cartesian (return z coordinate)
double GetUz(double ur, double ut, double sintheta, double costheta);

// Convert a vector from Cartesian to polar (return r coordinate)
double GetUr(double ux, double uy, double uz, double sintheta, double costheta,
   double sinphi, double cosphi);

// Convert a vector from Cartesian to polar (return theta coordinate)
double GetUt(double ux, double uy, double uz, double sintheta, double costheta,
   double sinphi, double cosphi);

// Convert a vector from Cartesian to polar (return phi coordinate)
double GetUp(double ux, double uy, double sinphi, double cosphi);

// Convert a vector from polar to Cartesian
void PolarToCart(const double *polar, double *cart, double sintheta,
   double costheta, double sinphi, double cosphi);

// Convert a vector from Cartesian to polar
void CartToPolar(const double *cart, double *polar, double sintheta,
   double costheta, double sinphi, double cosphi);

// Theta (polar) angle of a vector from the origin
double GetT(double *v);

// Theta (polar) angle of a vector from the origin - component input
double GetT(double nx, double ny, double nz);

// Phi (azimuthal) angle of a vector from the origin
double GetP(double *v);

// Phi (azimuthal) angle of a vector from the origin - component input
double GetP(double nx, double ny);


//==============================================================================
// Spherical trigonometry
//==============================================================================


// Cosine of the distance along a great circle connecting two points
double SphericalDist(double lat1, double long1, double lat2, double long2);

// Area of a spherical triangle on a unit sphere
double SphTriArea(double lat1, double long1, double lat2, double long2,
   double lat3, double long3);

// Calculate a point of intersection of two great circles - vector version
void GreatCircleInt(double *v1, double *v2, double *v3, double *v4, double *v5);

// Calculate a point of intersection of two great circles - lat/lon version
void GreatCircleInt(double lat1, double long1, double lat2, double long2,
   double lat3, double long3, double lat4, double long4, double &lat5,
   double &long5);

// Find the center of mass of a spherical triangle
void MassCenter(double *v1, double *v2, double *v3, double *cm);

// Compute a vector bisecting the angle between two vectors
void Bisect(double *v1, double *v2, double *v3);

// Compute a direction bisecting a given angle
void Bisect(double lat1, double long1, double lat2, double long2,
   double &lat3, double &long3);

// Polar angle of a point on a tilted great circle
double TiltedCircle(double alpha, double phi);

// Search an array for the point closest to a given point on a sphere
int NearestPoint(int npoints, double *plat, double *plong, double theta,
   double phi);

// Calculate x and y of a point in the local coordinate system
void PointCoordPlanar(double unit, double lat1, double long1, double lat2,
   double long2, double &x, double &y);


//==============================================================================
// Vector algebra
//==============================================================================


// Copy one vector to another
inline void Copy(const double *v1, double *v2) {memcpy(&v2[1], &v1[1], 3 * sizeof(double));};
   
// Compute the norm of a vector
double Norm(const double *v);

// Compute the norm of a complex vector
double NormComplex(const std::complex<double> *v);

// Compute the square of the norm
double Norm2(const double *v);

// Normalize a vector
void Normalize(double *v);

// Reflect a vector
void Reflect(double *v);

// Compute a scalar product
double ScalarProduct(const double *v1, const double *v2);

// Compute a vector product
void VectorProduct(const double *v1, const double *v2, double *v3);

// Compute a vector product of complex vectors
void VectorProductComplex(const std::complex<double> *v1, const std::complex<double> *v2, std::complex<double> *v3);

// Compute a triple product
double TripleProduct(const double *v1, const double *v2, const double *v3);

// Interpolate a variable using linear 3-point interpolation on a plane
double Interpolate3(double x1, double x2, double x3, double y1, double y2,
   double y3, double f1, double f2, double f3);

// Compute a normal vector to a plane defined with 3 unit vectors
void PlaneNormal(const double *v1, const double *v2, const double *v3,
   double *normal);


//==============================================================================
// Vector calculus
//==============================================================================


// Compute gradient of a scalar field
double *Gradient(int nn, double **normals, double area, double *elengths,
   double r, double *dr, double *var);

// Compute divergence of a vector field from the grid
double DivergenceC(int nn, double **normals, double area, double *elengths,
   double r, double *dr, double **var);

// Compute divergence of a vector field from face values
double DivergenceF(int nn, double **normals, double area, double *elengths,
   double r, double dr, double **var);

// Compute divergence of a vector field from normal face values
double DivergenceN(int nn, double area, double *elengths, double r, double dr,
   double *var);

// Compute curl of a vector field from the grid
double *Curl(int nn, double **normals, double area, double *elengths, double r,
   double *dr, double **var);

#endif


//==============================================================================
// End of the file "geo_coord.hh"
//==============================================================================
