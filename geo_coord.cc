//==============================================================================
// Start of the file "geo_coord.cc"
//
// This file is part of the SPECTRUM suite version 5a.
//
// Written by Vladimir Florinski
//==============================================================================


#include "geo_memory.hh"
#include "geo_coord.hh"

using namespace std;


//==============================================================================
// Cartesian-polar coordinate transformation
//==============================================================================


//------------------------------------------------------------------------------
// Convert a vector from polar to Cartesian (return x coordinate)
// Input:  ur - radial component
// Input:  ut - theta component
// Input:  up - phi component
// Input:  sintheta - sin(theta)
// Input:  costheta - cos(theta)
// Input:  sinphi - sin(phi)
// Input:  cosphi - cos(phi)
// Return: x-component
//------------------------------------------------------------------------------

double GetUx(double ur, double ut, double up, double sintheta, double costheta,
   double sinphi, double cosphi)
{
   return (ur * sintheta + ut * costheta) * cosphi - up * sinphi;
};


//------------------------------------------------------------------------------
// Convert a vector from polar to Cartesian (return y coordinate)
// Input:  ur - radial component
// Input:  ut - theta component
// Input:  up - phi component
// Input:  sintheta - sin(theta)
// Input:  costheta - cos(theta)
// Input:  sinphi - sin(phi)
// Input:  cosphi - cos(phi)
// Return: y-component
//------------------------------------------------------------------------------

double GetUy(double ur, double ut, double up, double sintheta, double costheta,
   double sinphi, double cosphi)
{
   return (ur * sintheta + ut * costheta) * sinphi + up * cosphi;
};


//------------------------------------------------------------------------------
// Convert a vector from polar to Cartesian (return z coordinate)
// Input:  ur - radial component
// Input:  ut - theta component
// Input:  sintheta - sin(theta)
// Input:  costheta - cos(theta)
// Return: z-component
//------------------------------------------------------------------------------

double GetUz(double ur, double ut, double sintheta, double costheta)
{
   return ur * costheta - ut * sintheta;
};


//------------------------------------------------------------------------------
// Convert a vector from Cartesian to polar (return r coordinate)
// Input:  ux - x component
// Input:  uy - y component
// Input:  uz - z component
// Input:  sintheta - sin(theta)
// Input:  costheta - cos(theta)
// Input:  sinphi - sin(phi)
// Input:  cosphi - cos(phi)
// Return: r-component
//------------------------------------------------------------------------------

double GetUr(double ux, double uy, double uz, double sintheta, double costheta,
   double sinphi, double cosphi)
{
   return (ux * cosphi + uy * sinphi) * sintheta + uz * costheta;
};


//------------------------------------------------------------------------------
// Convert a vector from Cartesian to polar (return theta coordinate)
// Input:  ux - x component
// Input:  uy - y component
// Input:  uz - z component
// Input:  sintheta - sin(theta)
// Input:  costheta - cos(theta)
// Input:  sinphi - sin(phi)
// Input:  cosphi - cos(phi)
// Return: theta-component
//------------------------------------------------------------------------------

double GetUt(double ux, double uy, double uz, double sintheta, double costheta,
   double sinphi, double cosphi)
{
   return (ux * cosphi + uy * sinphi) * costheta - uz * sintheta;
};


//------------------------------------------------------------------------------
// Convert a vector from Cartesian to polar (return phi coordinate)
// Input:  ux - x component
// Input:  uy - y component
// Input:  sinphi - sin(phi)
// Input:  cosphi - cos(phi)
// Return: phi-component
//------------------------------------------------------------------------------

double GetUp(double ux, double uy, double sinphi, double cosphi)
{
   return -ux * sinphi + uy * cosphi;
};


//------------------------------------------------------------------------------
// Convert a vector from polar to Cartesian
// Input:  polar[4] - vector in polar coordinates
// Output: cart[4] - vector in Cartesian coordinates
// Input:  sintheta - sin(theta)
// Input:  costheta - cos(theta)
// Input:  sinphi - sin(phi)
// Input:  cosphi - cos(phi)
//------------------------------------------------------------------------------

void PolarToCart(const double *polar, double *cart, double sintheta,
   double costheta, double sinphi, double cosphi)
{
   double urt = polar[1] * sintheta + polar[2] * costheta;
   cart[1] = urt * cosphi - polar[3] * sinphi;
   cart[2] = urt * sinphi + polar[3] * cosphi;
   cart[3] = polar[1] * costheta - polar[2] * sintheta;
};


//------------------------------------------------------------------------------
// Convert a vector from Cartesian to polar
// Intput: cart[4] - vector in Cartesian coordinates
// Ounput: polar[4] - vector in polar coordinates
// Input:  sintheta - sin(theta)
// Input:  costheta - cos(theta)
// Input:  sinphi - sin(phi)
// Input:  cosphi - cos(phi)
//------------------------------------------------------------------------------
void CartToPolar(const double *cart, double *polar, double sintheta,
   double costheta, double sinphi, double cosphi)
{
   double uxy = cart[1] * cosphi + cart[2] * sinphi;
   polar[1] = uxy * sintheta + cart[3] * costheta;
   polar[2] = uxy * costheta - cart[3] * sintheta;
   polar[3] = -cart[1] * sinphi + cart[2] * cosphi;
};


//------------------------------------------------------------------------------
// Theta (polar) angle of a vector from the origin
// Input:  v[4] - input vector
// Return: theta angle
//------------------------------------------------------------------------------

double GetT(double *v)
{
   double theta = atan2(sqrt(Sqr(v[1]) + Sqr(v[2])), v[3]);
   if(theta < 0.0) theta += twopi;
   return theta;
};


//------------------------------------------------------------------------------
// Theta (polar) angle of a vector from the origin - component input
// Input:  nx - x component
// Input:  ny - y component
// Input:  nz - z component
// Return: theta angle
//------------------------------------------------------------------------------

double GetT(double nx, double ny, double nz)
{
   double theta = atan2(sqrt(Sqr(nx) + Sqr(ny)), nz);
   if(theta < 0.0) theta += twopi;
   return theta;
};


//------------------------------------------------------------------------------
// Phi (azimuthal) angle of a vector from the origin
// Input:  v[4] - input vector
// Return: phi angle
//------------------------------------------------------------------------------

double GetP(double *v)
{
   double phi = atan2(v[2], v[1]);
   if(phi < 0.0) phi += twopi;
   return phi;
};


//------------------------------------------------------------------------------
// Phi (azimuthal) angle of a vector from the origin - component input
// Input:  nx - x component
// Input:  ny - y component
// Return: phi angle
//------------------------------------------------------------------------------

double GetP(double nx, double ny)
{
   double phi = atan2(ny, nx);
   if(phi < 0.0) phi += twopi;
   return phi;
};


//==============================================================================
// Spherical trigonometry
//==============================================================================


//------------------------------------------------------------------------------
// Cosine of the distance along a great circle connecting two points
// Input:  lat1 - latitude of the first point
// Input:  long1 - longitude of the first point
// Input:  lat2 - latitude of the second point
// Input:  long2 - longitude of the second point
// Return: cosine of the angle point 1 - origin - point 2
//------------------------------------------------------------------------------

double SphericalDist(double lat1, double long1, double lat2, double long2)
{
   return sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(long1 - long2);
};


//------------------------------------------------------------------------------
// Area of a spherical triangle on a unit sphere
// Input:  lat1 - latitude of the first point
// Input:  long1 - longitude of the first point
// Input:  lat2 - latitude of the second point
// Input:  long2 - longitude of the second point
// Input:  lat3 - latitude of the third point
// Input:  long3 - longitude of the third point
// Return: triangle area
//------------------------------------------------------------------------------

double SphTriArea(double lat1, double long1, double lat2, double long2,
   double lat3, double long3)
{
   double cosa, cosb, cosc, sina, sinb, sinc, A, B, C;

// Calculate the sides of the triangles
   cosa = SphericalDist(lat2, long2, lat3, long3);
   cosb = SphericalDist(lat1, long1, lat3, long3);
   cosc = SphericalDist(lat1, long1, lat2, long2);
   sina = sqrt(1.0 - cosa * cosa);
   sinb = sqrt(1.0 - cosb * cosb);
   sinc = sqrt(1.0 - cosc * cosc);

// Calculate the angles of the triangle using spherical law of cosines
   A = acos((cosa - cosb * cosc) / (sinb * sinc));
   B = acos((cosb - cosa * cosc) / (sina * sinc));
   C = acos((cosc - cosa * cosb) / (sina * sinb));

// Spherical excess
   return A + B + C - M_PI;
};


//------------------------------------------------------------------------------
// Calculate a point of intersection of two great circles. Formula may be found
// at http://geospatialmethods.org/spheres/GCIntersect.html#GCIGC
// Input:  v1 - unit vector to point 1 on the first circle
// Input:  v2 - unit vector to point 2 on the first circle
// Input:  v3 - unit vector to point 1 on the second circle
// Input:  v4 - unit vector to point 2 on the second circle
// Output: v5 - unit vector to the intersection point
//------------------------------------------------------------------------------

void GreatCircleInt(double *v1, double *v2, double *v3, double *v4, double *v5)
{
   double v1x2[4], v3x4[4], g, h, k;

// vector product r1 x r2 defines the first plane
   VectorProduct(v1, v2, v1x2);

// vector product r3 x r4 defines the second plane
   VectorProduct(v3, v4, v3x4);

// Find the line of intersection between the planes
   h =  (v3x4[1] * v1x2[3] - v3x4[3] * v1x2[1])
      / (v3x4[2] * v1x2[1] - v3x4[1] * v1x2[2]);
   g = -(v1x2[2] * h + v1x2[3]) / v1x2[1];
   k = 1.0 / sqrt(Sqr(g) + Sqr(h) + 1.0);

   v5[1] = g * k;
   v5[2] = h * k;
   v5[3] = k;

// Choose the direction that is closest to point 1 (angle < 90 deg).
   if(ScalarProduct(v1, v5) < 0.0) Reflect(v5);
   Normalize(v5);
};


//------------------------------------------------------------------------------
// Second version of the great circle intersection point finder using latitude-
// longitude notation.
// Input:  lat1 - latitude of point 1 on the first circle
// Input:  long1 - longitude of point 1 on the first circle
// Input:  lat2 - latitude of point 2 on the first circle
// Input:  long2 - longitude of point 2 on the first circle
// Input:  lat3 - latitude of point 1 on the second circle
// Input:  long3 - longitude of point 1 on the second circle
// Input:  lat4 - latitude of point 2 on the second circle
// Input:  long4 - longitude of point 3 on the second circle
// Output: lat5 - latitude of the intersection point
// Output: long5 - longitude of the intersection point
//------------------------------------------------------------------------------

void GreatCircleInt(double lat1, double long1, double lat2, double long2,
   double lat3, double long3, double lat4, double long4, double &lat5,
   double &long5)
{
   double v1[4], v2[4], v3[4], v4[4], v5[4];

// Cartesian coordinates of the points
   PolarToCart(radial_unit_vec, v1, cos(lat1), sin(lat1), sin(long1), cos(long1));
   PolarToCart(radial_unit_vec, v2, cos(lat2), sin(lat2), sin(long2), cos(long2));
   PolarToCart(radial_unit_vec, v3, cos(lat3), sin(lat3), sin(long3), cos(long3));
   PolarToCart(radial_unit_vec, v4, cos(lat4), sin(lat4), sin(long4), cos(long4));

// Use the vector version to find the intersection point
   GreatCircleInt(v1, v2, v3, v4, v5);

// Convert to latitude-longitude
   lat5 = pitwo - GetT(v5);
   long5 = GetP(v5);
};


//------------------------------------------------------------------------------
// Find the center of mass of a spherical triangle given by its Cartesian
// components
// Input:  v1[4] - first vertex
// Input:  v2[4] - second vertex
// Input:  v3[4] - third vertex
// Output: cm[4] - center of mass (on the sphere)
//------------------------------------------------------------------------------
void MassCenter(double *v1, double *v2, double *v3, double *cm)
{
   int i;
   double v1x2[4], v2x3[4], v3x1[4], Ah, Bh, Ch;

// Compute vectors normal to the sides of the triangle on the sphere
   VectorProduct(v1, v2, v1x2);
   VectorProduct(v2, v3, v2x3);
   VectorProduct(v3, v1, v3x1);
   Normalize(v1x2);
   Normalize(v2x3);
   Normalize(v3x1);

// Compute angles between the vectors
   Ah = ScalarProduct(v2, v3);
   Bh = ScalarProduct(v3, v1);
   Ch = ScalarProduct(v1, v2);

   for(i = 1; i <= 3; i++) cm[i] = v1x2[i] * Ch + v2x3[i] * Ah + v3x1[i] * Bh;
   Normalize(cm);
};


//------------------------------------------------------------------------------
// Compute a unit vector bisecting the angle between two given vectors
// (user is responsible for normalizing the input vectors)
// Input:  v1 - first vector
// Input:  v2 - second vector
// Output: v3 - bisecting vector
//------------------------------------------------------------------------------

void Bisect(double *v1, double *v2, double *v3)
{
   for(int i = 1; i <= 3; i++) v3[i] = (v1[i] + v2[i]) / 2.0;
   Normalize(v3);
};


//------------------------------------------------------------------------------
// Calculate a direction bisecting an angle between directions vec1 and vec2
// Input:  lat1 - latitude of point 1
// Input:  long1 - longitude of point 1
// Input:  lat2 - latitude of point 2
// Input:  long2 - longitude of point 2
// Output: lat3 - latitude of the bisection
// Output: long3 - longitude of the bisection
//------------------------------------------------------------------------------

void Bisect(double lat1, double long1, double lat2, double long2,
   double &lat3, double &long3)
{
   double v1[4], v2[4], bisection[4];

// Cartesian coordinates of the twp points
   v1[1] = cos(lat1) * cos(long1);
   v1[2] = cos(lat1) * sin(long1);
   v1[3] = sin(lat1);

   v2[1] = cos(lat2) * cos(long2);
   v2[2] = cos(lat2) * sin(long2);
   v2[3] = sin(lat2);

// find the midpoint and normalize it
   for(int i = 1; i <= 3; i++) bisection[i] = (v1[i] + v2[i]) / 2.0;
   Normalize(bisection);
   lat3 = pitwo - GetT(bisection);
   long3 = GetP(bisection);
};


//------------------------------------------------------------------------------
// Polar angle of a point on a great circle tilted by the angle alpha relative
// to the equatorial plane
// Input:  alpha - tilt angle
// Input:  phi - azimuthal angle
// Return: polar angle of the point
//------------------------------------------------------------------------------

double TiltedCircle(double alpha, double phi)
{
   return asin(sin(alpha) * sin(phi) / (fabs(cos(phi)) * sqrt(Sqr(cos(alpha))
      + Sqr(tan(phi)))));
};


//------------------------------------------------------------------------------
// Search an array for the point closest to a given point on a sphere
// Input:  npoints - size of the array
// Input:  plat[] - aray of latitudes
// Input:  plong[] - array of longitudes
// Input:  theta - theta of the point
// Input:  phi - phi of the point
// Return: index of the nearest point in the array
//------------------------------------------------------------------------------

int NearestPoint(int npoints, double *plat, double *plong, double theta,
   double phi)
{
   double cos_dist, cos_nearest, lat;
   int point, pnear = 1;
   lat = pitwo - theta;
   cos_nearest = tiny;

// Loop over all entries - inefficient for frequent use
   for(point = 1; point <= npoints; point++) {
      cos_dist = SphericalDist(plat[point], plong[point], lat, phi);
      if(cos_dist > cos_nearest) {
         cos_nearest = cos_dist; // nearest when cos angle is largest
         pnear = point;
      };
   };
   return pnear;
};


//------------------------------------------------------------------------------
// Calculate x and y of a point (lat2, long2) in the local coordinate system on
// the sphere centered at (lat1, long1). This routine is used to perform
// reconstruction on the sphere.
// Input:  unit - distance normalization constant (should be ~ 1/gridlevel)
// Input:  lat1 - latitude of the origin
// Input:  long1 - longitude of the origin
// Input:  lat2 - latitude of the point
// Input:  long2 - longitude of the point
// Output: x - local x coordinate of the point
// Output: y - local y coordinate of the point
//------------------------------------------------------------------------------

void PointCoordPlanar(double unit, double lat1, double long1, double lat2,
   double long2, double &x, double &y)
{
   double latC, longC, cosa, cosb, cosc, cosA, dist;

// Define point C just north of the origin to establish x-direction (N-S)
   latC = lat1 - unit; // unit must be small enough
   longC = long1;

// Fix if point is beyond the south pole
   if(latC < -pitwo) {
      latC = -M_PI - latC;
      longC += M_PI;
      if(longC > twopi) longC -= twopi;
   };

// Find the lengths of sides of the triangle between the origin, the point, and
// point C
   cosa = SphericalDist(latC, longC, lat2, long2);
   cosb = SphericalDist(latC, longC, lat1, long1);
   cosc = SphericalDist(lat1, long1, lat2, long2);
   dist = acos(cosc);

// Calculate the angles of the triangle
   if(1.0 - cosb < tiny || 1.0 - cosc < tiny) cosA = 1.0;
   else cosA = (cosa - cosb * cosc) / sqrt((1.0 - cosb * cosb)
                                         * (1.0 - cosc * cosc));
   if(cosA > 1.0) cosA = 1.0;
   if(cosA < -1.0) cosA = -1.0;

// Find coordinates as in planar geometry
   x = dist * cosA;
   y = dist * sqrt(1.0 - cosA * cosA);

// Check for sign
   if(long2 < long1) y = -y;
   if(long1 < M_PI && long2 > long1 + M_PI) y = -y;
   if(long1 > M_PI && long2 < long1 - M_PI) y = -y;
};


//==============================================================================
// Vector algebra
//==============================================================================


//------------------------------------------------------------------------------
// Compute the norm of a vector
// Input:  v[] - input vector
// Return: norm of v[]
//------------------------------------------------------------------------------

double Norm(const double *v)
{
   return sqrt(fmax(Sqr(v[1]) + Sqr(v[2]) + Sqr(v[3]), 0.0));
};


//------------------------------------------------------------------------------
// Compute the square of the norm
// Input:  v[] - input vector
// Return: norm² of v[]
//------------------------------------------------------------------------------

double Norm2(const double *v)
{
   return Sqr(v[1]) + Sqr(v[2]) + Sqr(v[3]);
};


//------------------------------------------------------------------------------
// Normalize a verctor
// Inout:  v[] - vector to normalize
//------------------------------------------------------------------------------

void Normalize(double *v)
{
   int i;
   double norm;
   norm = Norm(v);
   for(i = 1; i <= 3; i++) {
      v[i] /= norm;

// in case the initial vector is zero
      if(!isnormal(v[i])) v[i] = 0.0;
   };
};


//------------------------------------------------------------------------------
// Reflect a vector
// Inout:  v[] - vector to replace with its opposite
//------------------------------------------------------------------------------

void Reflect(double *v)
{
   for(int i = 1; i <= 3; i++) v[i] = -v[i];
};


//------------------------------------------------------------------------------
// Compute a scalar product of two vectors
// Input:  v1[] - first vector
// Input:  v2[] - second vector
// Return: scalar product 𝐯₁ ⋅ 𝐯₂
//------------------------------------------------------------------------------

double ScalarProduct(const double *v1, const double *v2)
{
   return v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
};


//------------------------------------------------------------------------------
// Compute a vector product of two vectors
// Input:  v1[] - first vector
// Input:  v2[] - second vector
// Output: v3[] - vector product 𝐯₁ ⨯ 𝐯₂
//------------------------------------------------------------------------------

void VectorProduct(const double *v1, const double *v2, double *v3)
{
   v3[1] = v1[2] * v2[3] - v1[3] * v2[2];
   v3[2] = v1[3] * v2[1] - v1[1] * v2[3];
   v3[3] = v1[1] * v2[2] - v1[2] * v2[1];
};


//------------------------------------------------------------------------------
// Compute a triple product
// Input:  v1[] - first vector
// Input:  v2[] - second vector
// Input:  v3[] - third vector
// Return: triple product 𝐯₁ ⋅ (𝐯₂ ⨯ 𝐯₃)
//------------------------------------------------------------------------------

double TripleProduct(const double *v1, const double *v2, const double *v3)
{
   double v2x3[4];
   VectorProduct(v2, v3, v2x3);
   return v1[1] * v2x3[1] + v1[2] * v2x3[2] + v1[3] * v2x3[3];
};


//------------------------------------------------------------------------------
// Interpolate a variable using linear 3-point interpolation on a plane
// Input:  x1 - x coordinate of point 1
// Input:  x2 - x coordinate of point 2
// Input:  x3 - x coordinate of point 3
// Input:  y1 - y coordinate of point 1
// Input:  y2 - y coordinate of point 2
// Input:  y3 - y coordinate of point 3
// Input:  f1 - value of the variable at point 1
// Input:  f2 - value of the variable at point 2
// Input:  f3 - value of the variable at point 3
// Return: interpolated value
//------------------------------------------------------------------------------

double Interpolate3(double x1, double x2, double x3, double y1, double y2,
   double y3, double f1, double f2, double f3)
{
   double a, b, c;
   a = x2 * y3 - x3 * y2;
   b = x3 * y1 - x1 * y3;
   c = x1 * y2 - x2 * y1;
   return (f1 * a + f2 * b + f3 * c) / (a + b + c);
};


//------------------------------------------------------------------------------
// Return a normal to a plane defined by three unit vectors
// Input:  v1[] - first vector
// Input:  v2[] - second vector
// Input:  v3[] - third vector
// Output: normal unit vector
//------------------------------------------------------------------------------

void PlaneNormal(const double *v1, const double *v2, const double *v3,
   double *normal)
{
   int i;
   double v2mv1[4], v3mv1[4];

   for(i = 1; i <= 3; i++) {
      v2mv1[i] = v2[i] - v1[i];
      v3mv1[i] = v3[i] - v1[i];
   };
   VectorProduct(v2mv1, v3mv1, normal);
   Normalize(normal);
   
// check for sign - normal must be in the same direction as the three vectors
   if(ScalarProduct(normal, v1) < 0.0) Reflect(normal);
};


//==============================================================================
// Vector calculus
//==============================================================================


//------------------------------------------------------------------------------
// Compute gradient of a scalar field inside a prism using Gauss theorem
// Input:  nn - number of sides of the polygon
// Input:  normals[4][nn+1] - array of *outward* normals (0th element is the
//         radial direction)
// Input:  area - area of the top/bottom surface on unit sphere
// Input:  elengths[nn+1] - lengths of prism edges on unit sphere
// Input:  r - radial distance to the center of the prism
// Input:  dr[4] - heights of the prisms l-1, l, l+1
// Input:  var[nn+3] - array of variables (0th element is the cell itself,
//         nn+1 and nn+2 are values below and above, respectively
// Return: gradient of scalar variable as a static array
//------------------------------------------------------------------------------

double *Gradient(int nn, double **normals, double area, double *elengths,
   double r, double *dr, double *var)
{
   int i, ie;
   static double gradu[4];
   double rdr, r2p, r2m;
   memset(gradu, 0, 4 * sizeof(double));
   rdr = r * dr[2];
   r2p = Sqr(r + dr[2] / 2.0);
   r2m = Sqr(r - dr[2] / 2.0);

// contributions from the top and bottom
   for(i = 1; i <= 3; i++) {
      gradu[i] += normals[i][0] * area * ((var[0] + (var[nn + 2] - var[0])
         * dr[2] / (dr[2] + dr[3])) * r2p - (var[0] + (var[nn + 1] - var[0])
         * dr[2] / (dr[2] + dr[1])) * r2m);
   };

// contributions from the sides
   for(ie = 1; ie <= nn; ie++) {
      for(i = 1; i <= 3; i++) {
         gradu[i] += normals[i][ie] * elengths[ie] * (var[0] + var[ie])
            * rdr / 2.0;
      };
   };
   
// divide by the cell volume
   for(i = 1; i <= 3; i++) gradu[i] /= (area * r * rdr);
   return gradu;
};


//------------------------------------------------------------------------------
// Compute divergence of a vector field inside a prism using Gauss theorem
// Input:  nn - number of sides of the polygon
// Input:  normals[4][nn+1] - array of *outward* normals (0th element is the
//         radial direction)
// Input:  area - spherical area of the top/bottom surface on unit sphere
// Input:  elengths[nn+1] - spherical lengths of prism edges on unit sphere
// Input:  r - radial distance of the center of the prism
// Input:  dr[4] - heights of the prisms l-1, l, l+1
// Input:  var[4][nn+3] - array of variables (0th element is the cell itself,
//         nn+1 and nn+2 are values below and above, respectively
// Return: divergence of vector variable
//------------------------------------------------------------------------------

double DivergenceC(int nn, double **normals, double area, double *elengths,
   double r, double *dr, double **var)
{
   int i, ie;
   double divu, rdr, r2p, r2m;
   divu = 0.0;
   rdr = r * dr[2];
   r2p = Sqr(r + dr[2] / 2.0);
   r2m = Sqr(r - dr[2] / 2.0);

// contributions from the top and bottom
   for(i = 1; i <= 3; i++) {
      divu += normals[i][0] * area * ((var[i][0] + (var[i][nn + 2] - var[i][0])
                                       * dr[2] / (dr[2] + dr[3])) * r2p
                                    - (var[i][0] + (var[i][nn + 1] - var[i][0])
                                       * dr[2] / (dr[2] + dr[1])) * r2m);
   };

// Contributions from the sides. Because edges intersect the lines connecting
// cell centers at midpoints, we can take arithmetic average of the two values.
   for(ie = 1; ie <= nn; ie++) {
      for(i = 1; i <= 3; i++) {
         divu += normals[i][ie] * elengths[ie] * (var[i][0] + var[i][ie])
            * rdr / 2.0;
      };
   };
   
// divide by the cell volume
   divu /= (area * r * rdr);
   return divu;
};


//------------------------------------------------------------------------------
// Another version of the divergence routine using face variables
// Input:  nn - number of sides of the polygon
// Input:  normals[4][nn+1] - array of *outward* normals (0th element is the
//         radial direction)
// Input:  area - spherical area of the top/bottom surface on unit sphere
// Input:  elengths[nn+1] - spherical lengths of prism edges on unit sphere
// Input:  r - radial distance of the center of the prism
// Input:  dr - height of the prism
// Input:  var_c[4][nn+3] - array of variables (0th element is unused,
//         1 - nn are side faces, nn+1 and nn+2 are bottom and top faces
// Return: divergence of vector variable
//------------------------------------------------------------------------------

double DivergenceF(int nn, double **normals, double area, double *elengths,
   double r, double dr, double **var)
{
   int i, ie;
   double divu, rdr, r2p, r2m;
   divu = 0.0;
   rdr = r * dr;
   r2p = Sqr(r + dr / 2.0);
   r2m = Sqr(r - dr / 2.0);

// contributions from the top and bottom
   for(i = 1; i <= 3; i++) {
      divu += normals[i][0] * area * (var[i][nn + 2] * r2p - var[i][nn + 1] * r2m);
   };

// contributions from the sides
   for(ie = 1; ie <= nn; ie++) {
      for(i = 1; i <= 3; i++) {
         divu += normals[i][ie] * elengths[ie] * var[i][ie] * rdr;
      };
   };
   
// divide by the cell volume
   divu /= (area * r * rdr);
   return divu;
};


//------------------------------------------------------------------------------
// Yet another version of the divergence routine using *normal* face variables
// Input:  nn - number of sides of the polygon
// Input:  area - spherical area of the top/bottom surface on unit sphere
// Input:  elengths[nn+1] - spherical lengths of prism edges on unit sphere
// Input:  r - radial distance of the center of the prism
// Input:  dr - height of the prism
// Input:  var_c[nn+3] - array of variables (0th element is unused,
//         1 - nn are side faces, nn+1 and nn+2 are bottom and top faces
// Return: divergence of vector variable
//------------------------------------------------------------------------------

double DivergenceN(int nn, double area, double *elengths, double r, double dr,
   double *var)
{
   int ie;
   double divu, rdr, r2p, r2m;
   rdr = r * dr;
   r2p = Sqr(r + dr / 2.0);
   r2m = Sqr(r - dr / 2.0);

// Contributions from the top and bottom. Note that the normal component on the
// bottom face is in the -R direction.
   divu = area * (var[nn + 2] * r2p + var[nn + 1] * r2m);

// contributions from the sides
   for(ie = 1; ie <= nn; ie++) {
      divu += elengths[ie] * var[ie] * rdr;
   };
   
// divide by the cell volume
   divu /= (area * r * rdr);
   return divu;
};


//------------------------------------------------------------------------------
// Compute curl of a vector field inside a prism using Gauss theorem
// Input:  nn - number of sides of the polygon
// Input:  normals[nn][3] - array of *outward* normals (0th element is the
//         radial direction)
// Input:  area - spherical area of the top/bottom surface on unit sphere
// Input:  elengths - spherical lengths of prism edges on unit sphere
// Input:  r - radial distance to the center of the prism
// Input:  dr - height of the prism
// Input:  var_c[nn+3][3] - array of variables (0th element is the cell itself,
//         nn+1 and nn+2 are values below and above, respectively
// Return: curl of vector variable as a static array
//------------------------------------------------------------------------------

double *Curl(int nn, double **normals, double area, double *elengths,
   double r, double *dr, double **var)
{
   const int curl1[4] = {0, 2, 3, 1};
   const int curl2[4] = {0, 3, 1, 2};
   int i, ie;
   static double curlu[4];
   double rdr, r2p, r2m;
   memset(curlu, 0, 4 * sizeof(double));
   rdr = r * dr[2];
   r2p = Sqr(r + dr[2] / 2.0);
   r2m = Sqr(r - dr[2] / 2.0);

// contributions from the top and bottom
   for(i = 1; i <= 3; i++) {
      curlu[i] += normals[curl1[i]][0] * area
               * ((var[curl2[i]][0] + (var[curl2[i]][nn + 2] - var[curl2[i]][0])
               * dr[2] / (dr[2] + dr[3])) * r2p
                - (var[curl2[i]][0] + (var[curl2[i]][nn + 1] - var[curl2[i]][0])
               * dr[2] / (dr[2] + dr[1])) * r2m)
                - normals[curl2[i]][0] * area
               * ((var[curl1[i]][0] + (var[curl1[i]][nn + 2] - var[curl1[i]][0])
               * dr[2] / (dr[2] + dr[3])) * r2p
                - (var[curl1[i]][0] + (var[curl1[i]][nn + 1] - var[curl1[i]][0])
               * dr[2] / (dr[2] + dr[1])) * r2m);
   };

// Contributions from the sides. Because edges intersect the lines connecting
// cell centers at midpoints, we can take arithmetic average of the two values.
   for(ie = 1; ie <= nn; ie++) {
      for(i = 1; i <= 3; i++) {
         curlu[i] += normals[curl1[i]][ie] * elengths[ie]
                  * (var[curl2[i]][0] + var[curl2[i]][ie]) * rdr / 2.0
                   - normals[curl2[i]][ie] * elengths[ie]
                  * (var[curl1[i]][0] + var[curl1[i]][ie]) * rdr / 2.0;
      };
   };

// divide by the cell volume
   for(i = 1; i <= 3; i++) curlu[i] /= (area * r * rdr);
   return curlu;
};


//==============================================================================
// End of the file "geo_coord.cc"
//==============================================================================
