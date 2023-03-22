/*
* Author: Tucker Quering
* Last revised: 2 December 2022
* Sets up uniform distribution of particles with weights
* according to a specified distribution function.
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include "pic_common.hh"

using namespace std;

// File names
const string params_dat = "params.dat";
const string fields_out = "fields.out";
const string particles_out = "prt00.out";

// Maximum velocity for the spherical quiet start
const double vmax = 3.5;

// Bit reversal function for quiet start courtesy of Kaijun Liu's code
double reverse(int num, int n){
  double power, rev;
  int inum, iquot, irem;
  
  rev = 0.;
  inum = num;
  power = 1.;
  
  do {
    iquot = inum/n;
    irem = inum - n*iquot;
    power /= n;
    rev += irem*power;
    inum = iquot;
  } while (inum > 0);
  
  return (rev);
}

// 2D efficient array memory management
// From Dr. Florinski's geo_memory.hh
// Copied rather than including file as it includes other definitions which overlap with pic_common.hh
template <class T> T **Create2D(int n, int m)
{
   T **array = new T *[n + 1];
   array[0] = new T[n * m + 1]; array[1] = array[0];
   for(int i = 2; i <= n; i++) array[i] = array[i - 1] + m;
   return array;
};
template <class T> void Delete2D(T **array) {delete[] array[0]; delete[] array;};

// Computes the value of a Maxwellian distribution
double compute_maxwellian(double density, double vt, double v_magnitude) {
    return (density / (Cube(sqrtpi*vt))) * exp(-Sqr(v_magnitude)/Sqr(vt));
}

// Structure holding functions / values for computing erf / inverse error function
// Courtesy of Numerical Recipes 3rd ed., (p. 264 - 265)
struct Erf {
    static const int ncof = 28;
    static const double cof[28];

    // For computing the error function
    inline double erf(double x) {
        if (x >= 0.) return 1.0 - erfccheb(x);
        else return erfccheb(-x) -1.0;
    }

    // For computing the complementary error function
    // This one can be called outside the struct
    inline double erfc(double x) {
        if (x >= 0.) return erfccheb(x);
        else return erfccheb(-x) - 1.0;
    }

    // Computes complementary erf using Chebyshev coefficients 
    // Should be called by erf() or erfc(), not by user
    double erfccheb(double z) {
        int j;
        double t, ty, tmp, d=0., dd=0.;
        if (z < 0.) throw("erfccheb requirs nonnegative argument");
        t = 2./(2. + z);
        ty = 4.*t - 2.;
        for (j = ncof-1; j > 0; j--) {
            tmp = d;
            d = ty*d - dd + cof[j];
            dd = tmp;
        }
        return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
}


    // Computes the inverse of the complementary error function
    double inverfc(double p) {
        double x, err, t, pp;
        if (p >= 2.0) return -100.; // Handles excessively large or small inputs
        if (p < 0.0) return 100.;
        pp = (p < 1.0) ? p : 2. - p;
        t = sqrt(-2.*log(pp/2.));
        x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
        for(int j = 0; j < 2; j++) {
            err = erfc(x) - pp;
            x += err/(1.12837916709551257*exp(-Sqr(x)) - x*err);
        }

        return (p < 1.0 ? x : -x);
    }

    // Computes the inverse error function
    double inverf(double p) {
        return inverfc(1. - p);
    }

};

const double Erf::cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
    1.9476473204185836e-2, -9.561514786808631e-3, -9.46595344482036e-4, 
    3.6683949785261e-4, 4.2523324806907e-5, -2.0278758112534e-5, 
    -1.624290004647e-6, 1.303655835580e-6, 1.5626441722e-8, -8.5328095915e-8,
    6.529054439e-9, 5.059343495e-9, -9.91364156e-10, -2.27365122e-10,
    9.6467911e-11, 2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13,
    -1.12708e-13, 3.81e-16, 7.106e-15, -1.523e-15, -9.4e-17, 1.21e-16, -2.8e-17};

// Struct for evaluating the integrated maxwellian we invert with the function below
// x0 = the value we would like to invert
// w0 = the value we evaluate the function we are finding roots of at.
struct mxw_int {
    double operator()(const double x0, const double w0) {
        return erf(w0) - (2*w0/sqrtpi)*exp(-Sqr(w0)) - x0;
    }

    double df(double w0) {
        return (4*Sqr(w0)/sqrtpi)*exp(-Sqr(w0));
    }
};

// Finds the root of the equation erf(w) - (2w/sqrt(pi))exp(-w^2) = x0 (detailed in notes)
// Uses a mis of bisection & Newton's method to ensure convergence - detaile on p. 460 in Numerical Recipes
double mxw_integral_inverse(double x0) {
    const int max_iterations = 100;
    Erf errorfunc = Erf();
    mxw_int funcd = mxw_int();
    double accuracy = 1.e-10;
    double xl = 5.0e-5;     // I hope to have complete convergence within this range of possible outputs
    double xh = vmax;       // If a value of x0 too low or high is passed to the function, we'll just return 0 or vmax respectively
    double fl = funcd(x0, xl);
    double fh = funcd(x0, xh);
    if (fl > 0.0 && fh > 0.0) return 0;     // Outlier cases; hopefully we won't see to many of these.
    if (fl < 0.0 && fh < 0.0) return 3.5;   // In theory there should be none on the upper bound
    // Worst-case scenario, if this happens too often, we'll simply re-introduce a singularity at the center
    // That's really not the end of the world in my book.

    double w0 = errorfunc.inverf(x0);
    double dxold = (xh - xl);
    double dx = dxold;
    double f = funcd(x0, w0);
    double df = funcd.df(w0);

    for(int j = 0; j < max_iterations; j++) {
        if ((((w0-xh)*df-f)*((w0-xl)*df-f) > 0.0) || (abs(2.0*f) > abs(dxold*df))) { // If Newton is out of range or failing to converge; bisect
            dxold = dx;
            dx = 0.5*(xh-xl);
            w0 = xl + dx;
            if (xl == w0) return w0;
        } else { // Otherwise, take a Newton step
            dxold = dx;
            dx = f/df;
            double temp = w0;
            w0 -= dx;
            if (temp == w0) return w0;
        }

        if (abs(dx) < accuracy) return w0; // Check if we have reached a suitably accurate value
        double f = funcd(x0, w0);
        double df = funcd.df(w0);
        if (f < 0.0) xl = w0;
        else xh = w0;
    }
    
    throw("mxw_integral_inverse hit the maximum number of iterations.");
}

int main(int argc, char *argv[]) {
    ifstream params_file;
    ofstream output_file;

    // Values to read from params file
    int Nspecies, Ncells, ppc;
    double xmax, B0, cos_theta, e_temp, density;

    // Directly computed values
    double v_t, wpi, wci, wpiwci, lambda_i, dx, Bx, By, zero;
    int field_cells, Nparts;

    // Values used in generating particles
    double **rvw;
    double w_sum, v, x, mu, theta, phi;

    string temps1, temps2;
    
    params_file.open(params_dat.c_str(), ifstream::in);

    params_file >> temps1 >> Nspecies;
    if (Nspecies > 1) {
        cerr << "WARNING: Only generating particles for species #1.\n.";
    }

    params_file >> temps1 >> Ncells;
    params_file >> temps1 >> xmax;
    params_file >> temps1 >> temps2; // dt is discarded as it is not needed

    params_file >> temps1 >> B0;
    params_file >> temps1 >> cos_theta;
    params_file >> temps1 >> e_temp;

    params_file >> temps1 >> temps2; // discard driving current data + distribution type
    params_file >> temps1 >> temps2;
    params_file >> temps1 >> temps2;

    params_file >> temps1 >> ppc;
    params_file >> temps1 >> density;

    v_t = sqrt(2 * kboltz * e_temp / p_mass);
    wpi = 2.0 * sqrtpi * echarge * sqrt(density / p_mass);
    wci = echarge * B0 / (p_mass * splight);
    wpiwci = wpi/wci;
    lambda_i = splight / wpi;
    dx = xmax / Ncells;

    Nparts = Ncells * ppc;

    cerr << hline << endl;
    cerr << "Computing for " << Ncells << " cells, with " << ppc << " particles per cell.\n";
    cerr << "Ion inertial length is " << lambda_i << " cm.\n";
    cerr << "Thermal speed is " << v_t << " cm/s.\n";
    cerr << "Plasma to cyclotron frequency ratio is " << wpiwci << ".\n";
    cerr << hline << endl;

    Bx = cos_theta / wpiwci;
    By = sqrt(1 - Sqr(cos_theta)) / wpiwci;
    zero = 0.0;
    field_cells = Ncells + 2; // Ghost cells

    // Write the field values first since these are simple and uniform
    output_file.open(fields_out.c_str(), ofstream::out | ofstream::binary);
    output_file.write((char *)&field_cells, sizeof(int));
    output_file.write((char *)&zero, sizeof(double)); // Corresponds to simtime in the file
    output_file.write((char *)&zero, sizeof(double)); // This one is one of the extraneous zeros in the file that I just can't explain.
    for(int n = 0; n < field_cells; n++) {
        output_file.write((char *)&Bx, sizeof(double));
        output_file.write((char *)&By, sizeof(double));
        output_file.write((char *)&zero, sizeof(double));
    }

    // Electric field is set to zero everywhere, plus the "extraneous zero" at the beginning.
    for(int n = 0; n < 3*field_cells + 1; n++) {
        output_file.write((char *)&zero, sizeof(double));
    }

    output_file.close();
    cerr << "Finished writing field data.\n";
    cerr << "Proceeding to generate " << Nparts << " particles.\n";

    rvw = Create2D<double>(Nparts, 7);

    for(int n = 1; n <= Nparts; n++) {
        if(!(n % (Nparts/10)) && n>0) {
            cerr << "Particle generation is " << (n*100)/Nparts << "\% complete.\n";
        }
        // Begin by setting up particle values in physical units
        // Uses a quiet start algorithm for a nice distribution in x, vx, vy, and vz
        // This setup ensures each cell receives the same number of particles
        rvw[n][1] = (floor((n-1)/ppc) + reverse(n-1, 2)) * dx * lambda_i; // x
        rvw[n][2] = 0.0; // y
        rvw[n][3] = 0.0; //z

        // Place velocities using a Maxwellian quiet start sphere
        x = reverse(n, 3) * 0.9999;
        mu = reverse(n, 5)*2 - 1;
        phi = reverse(n, 7)*twopi;

        // Convert the uniformly chosen x to a maxwellian velocity using the inverse function
        v = mxw_integral_inverse(x) * v_t;
        // And of course mu = cos(theta), so
        theta = acos(mu);

        // Switch from spherical coordinates to cartesian when writing to the array
        rvw[n][4] = v * sin(theta) * cos(phi);
        rvw[n][5] = v * sin(theta) * sin(phi);
        rvw[n][6] = v * cos(theta);

        // IN THEORY due to the way we have initialized the particles, we can return to uniform weights here.
        rvw[n][7] = 1; // 1 as a placeholder but hey, maybe we could just do 1/ppc and stop normalizing? That probably doesn't generalize to putting an arbitrary distribution here though.

        // Now that the weight has been computed, we switch from physical units to code units
        // Distance is measured in inertial lengths, velocity is measured relative to the speed of light
        rvw[n][1] /= lambda_i;
        rvw[n][4] /= splight;
        rvw[n][5] /= splight;
        rvw[n][6] /= splight;
    }

    // At this point, rvw should be fully populated with values
    // The weights of the particles must now be normalized so that the sum of the weights = Ncells
    w_sum = 0.0;
    for(int n = 1; n <= Nparts; n++) {
        w_sum += rvw[n][7];
    }

    for(int n = 1; n <= Nparts; n++) {
        rvw[n][7] *= ((double)Ncells/w_sum);
    }
    cerr << "Finished normalizing particle weights.\n";
    
    output_file.open(particles_out.c_str(), ofstream::out | ofstream::binary);
    output_file.write((char *)&Nparts, sizeof(int));
    output_file.write((char *)&rvw[0][1], (Nparts * 7) * sizeof(double));
    output_file.close();
    cerr << "Finished writing particle data.\n";
}