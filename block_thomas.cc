//==============================================================================
// Start of the file "block_thomas.cc"
//
// Version 3d.
//
// Written by Vladimir Florinski
//==============================================================================


#include "block_thomas.hh"

using namespace std;


double ZeroMatr[4] = {0.0, 0.0, 0.0, 0.0};
double UnitMatr[4] = {1.0, 0.0, 0.0, 1.0};


// Subtract two 2x1 vectors: V3=V1-V2
inline void Sub_VV2(double *V1, double *V2, double *V3)
{
   V3[0] = V1[0] - V2[0];
   V3[1] = V1[1] - V2[1];
};


// Add two 2x2 matrices: M3=M1+M2
inline void Add_MM2(double *M1, double *M2, double *M3)
{
   M3[0] = M1[0] + M2[0];
   M3[1] = M1[1] + M2[1];
   M3[2] = M1[2] + M2[2];
   M3[3] = M1[3] + M2[3];
};


// Subtract two 2x2 matrices M3=M1-M2
inline void Sub_MM2(double *M1, double *M2, double *M3)
{
   M3[0] = M1[0] - M2[0];
   M3[1] = M1[1] - M2[1];
   M3[2] = M1[2] - M2[2];
   M3[3] = M1[3] - M2[3];
};


// Multiply a 2x2 matrix by a 2x1 vector: V3=M1*V2
inline void Mul_MV2(double *M1, double *V2, double *V3)
{
   double V4[2];

// V3 could be the same as V2, so we need to multiply out of place
   V4[0] = M1[0] * V2[0] + M1[1] * V2[1];
   V4[1] = M1[2] * V2[0] + M1[3] * V2[1];
   V3[0] = V4[0];
   V3[1] = V4[1];
};


// Multiply two 2x2 matrices M3=M1*M2
inline void Mul_MM2(double *M1, double *M2, double *M3)
{
   double M4[4];

// M3 could be one of M1 or M2, so we need to multiply out of place
   M4[0] = M1[0] * M2[0] + M1[1] * M2[2];
   M4[1] = M1[0] * M2[1] + M1[1] * M2[3];
   M4[2] = M1[2] * M2[0] + M1[3] * M2[2];
   M4[3] = M1[2] * M2[1] + M1[3] * M2[3];

// this is probably faster than memcpy()
   M3[0] = M4[0];
   M3[1] = M4[1];
   M3[2] = M4[2];
   M3[3] = M4[3];
};


// Invert a 2x2 matrix: M2=M1^-1
inline void MInv_M2(double *M1, double *M2)
{
   double det = M1[0] * M1[3] - M1[1] * M1[2];
   M2[0] =  M1[3] / det;
   M2[1] = -M1[1] / det;
   M2[2] = -M1[2] / det;
   M2[3] =  M1[0] / det;
};


//------------------------------------------------------------------------------
// Three-diagonal periodic Thomas algorithm for block size 2
//------------------------------------------------------------------------------
// Input  | N        | size of the matrix
// Input  | a[][4]   | array of "a" (left of diagonal) matrices
// Input  | b[][4]   | array of "b" (diagonal) matrices
// Input  | c[][4]   | array of "c" (right of diagonal) matrices
// Input  | d[][2]   | right hand side vector
// Output | x[][2]   | solution vector - cannot point to "d"
//------------------------------------------------------------------------------

void ThomasPeriodicBlock2(int N, double **a, double **b, double **c, double **d,
                          double **x)
{
   int i;
   double binv[4], aa[4], ab[4], ac[4], bc[4], ca[4], cc[4], ad[2], cd[2],
          ax[2], cx[2];
   double bcinv[4], cbcinv[4], cbcinva[4], abcbcinva[4], abcbcinvainv[4],
          cbcinvd[2], dcbcinvd[2], dax[2];

// Normalize the first row
   MInv_M2(b[0], binv);
   Mul_MM2(binv, a[0], a[0]);
   Mul_MM2(binv, c[0], c[0]);
   Mul_MV2(binv, d[0], d[0]);

// Eliminate the lower diagonal with normalization
   for(i = 1; i < N; i++) {

// Compute b_i
      Mul_MM2(a[i], c[i - 1], ac);
      Sub_MM2(b[i], ac, b[i]);
      
// Compute d_i
      Mul_MV2(a[i], d[i - 1], ad);
      Sub_VV2(d[i], ad, d[i]);
      
// Compute a_i - last, because its current values are needed above
      Mul_MM2(a[i], a[i - 1], aa);
      Sub_MM2(ZeroMatr, aa, a[i]);

// Normalize the row. Note that we do not compute b_i explicitly - we assume it
// is a unit matrix from this point on.
      MInv_M2(b[i], binv);
      Mul_MM2(binv, a[i], a[i]);
      Mul_MM2(binv, c[i], c[i]);
      Mul_MV2(binv, d[i], d[i]);
   };

// Eliminate the upper diagonal
   for(i = N - 2; i >= 0; i--) {
      
// Compute a_i
      Mul_MM2(c[i], a[i + 1], ca);
      Sub_MM2(a[i], ca, a[i]);

// Compute d_i
      Mul_MV2(c[i], d[i + 1], cd);
      Sub_VV2(d[i], cd, d[i]);

// Compute c_i - last, because its current values are needed above
      Mul_MM2(c[i], c[i + 1], cc);
      Sub_MM2(ZeroMatr, cc, c[i]);
   };
      
// Find  x_1 and x_N
   Add_MM2(a[N - 1], UnitMatr, ab);
   Add_MM2(UnitMatr, c[0], bc);
   MInv_M2(bc, bcinv);
   Mul_MM2(c[N - 1], bcinv, cbcinv);
   Mul_MM2(cbcinv, a[0], cbcinva);
   Mul_MV2(cbcinv, d[0], cbcinvd);
   Sub_MM2(ab, cbcinva, abcbcinva);
   MInv_M2(abcbcinva, abcbcinvainv);
   Sub_VV2(d[N - 1], cbcinvd, dcbcinvd);
   Mul_MV2(abcbcinvainv, dcbcinvd, x[N - 1]);
   Mul_MV2(a[0], x[N - 1], ax);
   Sub_VV2(d[0], ax, dax);
   Mul_MV2(bcinv, dax, x[0]);
   
// Back-substitution, using the result that b_i are unit matrices
   for(i = 1; i < N - 1; i++) {
      Mul_MV2(a[i], x[N - 1], ax);
      Mul_MV2(c[i], x[0], cx);
      Sub_VV2(d[i], ax, x[i]);
      Sub_VV2(x[i], cx, x[i]);
   };
};


//==============================================================================
// End of the file "block_thomas.cc"
//==============================================================================
