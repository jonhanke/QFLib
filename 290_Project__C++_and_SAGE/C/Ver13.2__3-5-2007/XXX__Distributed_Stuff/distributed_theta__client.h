


////////////////////////////////////////////////////////////////
// Writes the boolean ternary theta function as a binary file //
// NOTE: This is just a copy of the routine in maketheta.cc!  //
// (with some minor modifications to the directory and filename...)
////////////////////////////////////////////////////////////////

void WriteTernaryThetaBinary_client(unsigned long theta_bin[], long form[6], unsigned long precision, char* hostname);




//////////////////////////////////////////////////////
// Find all vectors with Q(x) <= C in a given slice //
//////////////////////////////////////////////////////

void FastBinaryThetaDistributed(unsigned long theta[], double QQ[], double C, long slice_num); 



// Main program (server) for computing the distributed theta function

void MakeTernaryThetaDistributed__client(long QQ[6], unsigned long Ternary_Precision, long slice_number);



// -------------------------------------------

#if !defined(DISTRIBUTED_THETA__CLIENT_H)
#define DISTRIBUTED_THETA__CLIENT_H

#include "distributed_theta__client.cc"

#endif
