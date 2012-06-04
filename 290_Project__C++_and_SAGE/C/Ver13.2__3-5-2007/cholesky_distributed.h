

// Cholesky decomposition of a ternary form (both given by the 6 upper triangular elements).
// (NOTE: This was copied directly from maketheta.cc!)

void CholeskyTernary(long QQ[], double Cholesky[]);

void CholeskyDecomposition(Matrix_mpz QQ, double Cholesky[]);




// -------------------------------------------------



#if !defined(CHOLESKY_DISTRIBUTED_H)
#define CHOLESKY_DISTRIBUTED_H


#include "Matrix_mpz/Matrix_mpz.h"


#include "cholesky_distributed.cc"



/*
void __dummy_cholesky_distributed() {

  long QQ[1];
  double Cholesky[1];
  Matrix_mpz Q;

  CholeskyTernary(QQ, Cholesky);
  CholeskyDecomposition(Q, Cholesky);

}
*/


#endif
