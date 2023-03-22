/*!
\file pic_common.hh
\brief Declares some simple common constants and routines
\author Vladimir Florinski
*/

#ifndef PIC_COMMON
#define PIC_COMMON

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

//! Number of particle species
#define MAX_SPECIES 1

#define is_odd(i)   ((i) & 1)
#define is_even(i) !((i) & 1)

//! Number of ghost cells
#define GHOST 1

const std::complex<double> i(0.0, 1.0);

//! \f$2\pi\f$
const double twopi = 2.0 * M_PI;

//! \f$4\pi\f$
const double fourpi = 4.0 * M_PI;

//! \f$\pi/3\f$
const double piover3 = M_PI / 3.0;

//! \f$\sqrt{\pi}\f$
const double sqrtpi = sqrt(M_PI);

//! speed of light (cm/s)
const double splight = 2.99792458E+10;

//! elementary charge (esu)
const double echarge = 4.8032044E-10;  

//! mass of a proton (g)
const double p_mass = 1.6726219E-24;

//! mass of an electron (g)
const double e_mass = 9.10938356E-28;

// Boltzmann constant (erg/K)
const double kboltz = 1.3806505E-16;

//! 1 eV (erg)
const double one_ev = 1.60218E-12;

//! Polytropic index
const double pgamma = 5.0 / 3.0;

//! A horizontal line for output to the terminal
const std::string hline = "--------------------------------------------------------------------------------\n";

/*!
\brief Computes the square
\author Vladimir Florinski
\param[in] x The argument
\return \f$x^2\f$
*/
template <class T> inline T Sqr(T x) {return x * x;};

/*!
\brief Computes the cube
\author Vladimir Florinski
\param[in] x The argument
\return \f$x^3\f$
*/
template <class T> inline T Cube(T x) {return x * x * x;};

/*!
\brief Power of two
\author Vladimir Florinski
\param[in] n The exponent
\return \f$2^n\f$
*/
inline int Pow2(int n) {return 1 << n;};

/*!
\brief Scalar product
\author Vladimir Florinski
\param[in] v1 First vector
\param[in] v2 Second vector
\return \f$\mathbf{v}_1\cdot\mathbf{v}_2\f$
*/
template <class T> inline T ScalarProd(T* v1, T* v2)
{
   return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
};

/*!
\brief Magnitude of a vector
\author Vladimir Florinski
\param[in] v Vector
\return \f$|\mathbf{v}|\f$
*/
template <class T> inline T VectorMag(T* v)
{
   return sqrt(ScalarProd(v, v));
};

/*!
\brief Vector product
\author Vladimir Florinski
\param[in]  v1 First vector
\param[in]  v2 Second vector
\param[out] v3 \f$\mathbf{v}_1\times\mathbf{v}_2\f$
*/
template <class T> inline void VectorProd(T* v1, T* v2, T* v3)
{
   v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
   v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
   v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
};

/*!
\brief Normalize a vector
\author Vladimir Florinski
\param[in,out] v Vector
*/
template <class T> inline void Normalize(T* v)
{
   T norm = sqrt(ScalarProd(v, v));
   for(int xyz = 0; xyz < 3; xyz++) v[xyz] /= norm;
};

/*!
\brief Effective temperature based on characteristic velocity
\author Vladimir Florinski
\param[in] m Mass
\param[in] v Velocity
\return Effective temperature
*/
inline double EffectiveTemperature(double m, double v) {return m * Sqr(v) / (2.0 * kboltz);};

/*!
\brief Thermal speed based on temperature
\author Vladimir Florinski
\param[in] m Mass
\param[in] T Temperature
\return Thermal speed
*/
inline double ThermalSpeed(double m, double T) {return sqrt(2.0 * kboltz * T / m);};

/*!
\brief Plasma frequency
\author Vladimir Florinski
\param[in] m Mass
\param[in] q Charge
\param[in] n Number density
\return Plasma frequency
*/
inline double PlasmaFrequency(double m, double q, double n) {return fabs(q) * sqrt(fourpi * n / m);};

/*!
\brief Debye length
\author Vladimir Florinski
\param[in] q Charge
\param[in] n Number density
\param[in] n Temperature
\return Debye length
*/
inline double DebyeLength(double q, double n, double T) {return sqrt(kboltz * T / fourpi / n) / fabs(q);};

/*!
\brief Cyclotron frequency (non-relativistic)
\author Vladimir Florinski
\param[in] m Mass
\param[in] q Charge
\param[in] B magnetic field
\return Cyclotron frequency
*/
inline double CyclotronFrequency(double m, double q, double B) {return q * B / m / splight;};

/*!
\brief Larmor radius (non-relativistic)
\author Vladimir Florinski
\param[in] m Mass
\param[in] q Charge
\param[in] v Velocity
\param[in] B magnetic field
\return Larmor radius
*/
inline double LarmorRadius(double m, double q, double v, double B) {return m * v * splight / q / B;};

/*!
\brief Adiabatic pressure law
\author Vladimir Florinski
\param[in]  n0  Reference number density
\param[in]  n   Number density
\param[in]  T0  Reference temperature
\return Pressure
*/
inline double Pressure(double n0, double n, double T0) {return n0 * kboltz * T0 * pow(n / n0, pgamma);};

//! Color use for error messages
const std::string msg_color = "\033[31m";

//! Standard terminal text color
const std::string std_color = "\033[0m";

/*!
\brief Print an error message in color
\author Vladimir Florinski
\param[in] filename Source file name
\param[in] line     Source line number
\param[in] message  Message to print
*/
inline void PrintError(const char* filename, int line, const std::string& message)
{
   std::cerr << msg_color;
   std::cerr << filename << ":" << line << ": error: " << message << std::endl;
   std::cerr << std_color;
};

#endif
