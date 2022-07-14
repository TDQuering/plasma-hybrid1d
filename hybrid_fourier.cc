#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <complex.h>
#include <fftw3.h>

using namespace std;

const double xmax = 512.0;
const int Nx = 1024;
const double dx = xmax / Nx;

//const int Nt = 157; // 20 0 20
const int Nt = 785; // 200 100 200
const double dt = 0.02 * 40;


int main(void)
{
   int i, j, x, t, k;
   double kmax, omax, absv;
   double *spactime, *slice;
   fftw_complex *wavefreq, *spect;
   fftw_plan p1, p2;
   ifstream st_file;
   ofstream br_file, wf_file;

   kmax = M_PI / dx;
   omax = M_PI / dt;

   spactime = (double *)fftw_malloc(Nt * Nx * sizeof(double));
   slice = (double *)fftw_malloc(Nt * sizeof(double));
   
   wavefreq = (fftw_complex *)fftw_malloc(Nt * (Nx / 2 + 1) * sizeof(fftw_complex));
   spect = (fftw_complex *)fftw_malloc((Nt / 2 + 1) * sizeof(fftw_complex));

   p2 = fftw_plan_dft_r2c_2d(Nt, Nx, spactime, wavefreq, FFTW_ESTIMATE);
   p1 = fftw_plan_dft_r2c_1d(Nt, slice, spect, FFTW_ESTIMATE);

   st_file.open("spacetime.dat", ifstream::in | ifstream::binary);
   st_file.read((char *)spactime, Nt * Nx * sizeof(double));
   st_file.close();

   j = 500; // must be between 0 and Nx
   for(i = 0; i < Nt; i++) {
      k = i * Nx + j;
      slice[i] = spactime[k];
//      cout << i << " " << slice[i] << endl;
   };


// test: sine wave
/*
   for(i = 0; i < Nt; i++) {
      t = i * dt;
      for(j = 0; j < Nx; j++) {
         x = j * dx;
         k = i * Nx + j;
         spactime[k] = sin(1.0 * t);
      };
   };
*/

   fftw_execute(p2);
   fftw_execute(p1);
   fftw_destroy_plan(p2);
   fftw_destroy_plan(p1);

   for(i = 0; i < Nt / 2; i++) {
      cout << i * 2.0 * M_PI / (Nt * dt) << " " << cabs(spect[i]) << endl;
   };


   br_file.open("wavefreq.bov", ofstream::out);
   br_file << "DATA_FILE: wavefreq.dat\n";
   br_file << "DATA_SIZE: " << Nx / 4 << " " << Nt / 2 << " 1\n";
   br_file << "DATA_FORMAT: DOUBLE\n";
   br_file << "VARIABLE: AbsLambda\n";
   br_file << "DATA_ENDIAN: LITTLE\n";
   br_file << "CENTERING: zonal\n";
   br_file << "BRICK_ORIGIN: 0.0 0.0 0.0\n";
   br_file << "BRICK_SIZE: " << kmax / 2.0 << " " << omax << " 1.0\n";
   br_file.close();

   wf_file.open("wavefreq.dat", ofstream::out | ofstream::binary);
   for(i = 0; i < Nt / 2; i++) {

// cut off the large k part since this is strongly affected by the grid
      for(j = 0; j < Nx / 4; j++) {
         k = i * (Nx / 2 + 1) + j;
         absv = cabs(wavefreq[k]);
         wf_file.write((char *)&absv, sizeof(double));
      };
   };
   wf_file.close();

   fftw_free(spactime);
   fftw_free(slice);
   fftw_free(wavefreq);
   return 0;
};
