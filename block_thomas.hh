//==============================================================================
// Start of the file "block_thomas.hh"
//
// Version 3d.
//
// Written by Vladimir Florinski
//==============================================================================


#ifndef BLOCK_THOMAS
#define BLOCK_THOMAS


// Three-diagonal periodic Thomas algorithm for block size 2
void ThomasPeriodicBlock2(int N, double **a, double **b, double **c, double **d,
                          double **x);

#endif


//==============================================================================
// End of the file "block_thomas.hh"
//==============================================================================
