


using namespace std;

#include <iostream>


// Include the Matrix_mpz class
#include "Matrix_mpz/Matrix_mpz.h"




// Cholesky decomposition of a ternary form (both given by the 6 upper triangular elements).
// (NOTE: This was copied directly from maketheta.cc!)

void CholeskyTernary(long QQ[], double Cholesky[]) {

  const long n = 3;  // for ternary matrices

  double Q[n+1][n+1];  // Only the positive indices (1--3) are used... =)

  // 1. Initialize (from a symmetric matrix QQ) -- using 2 * QQ right now...
  long counter = 0;
  for (long i=1; i<=n; i++)
    for (long j=i; j<=n; j++) {
      if (i==j)
	Q[i][j] = 1 * QQ[counter];
      else 
	Q[i][j] = 0.5 * QQ[counter];
      counter++;
    }
  long i = 0;


  // 1a. Print the resulting form
  /*
    cout << endl;
    cout << "1a.  This is the starting (upper triangular) matrix of Q:" << endl;
    cout << "[ " << Q[1][1] << " " << Q[1][2] << " " << Q[1][3] << " ]" << endl;
    cout << "[ " << Q[2][1] << " " << Q[2][2] << " " << Q[2][3] << " ]" << endl;
    cout << "[ " << Q[3][1] << " " << Q[3][2] << " " << Q[3][3] << " ]" << endl;
    cout << endl;
  */
  
  
  // 2. Loop on i
  i = i + 1;
  while (i != n) {
    for (long j=i+1; j<=n; j++) {    
      Q[j][i] = Q[i][j];
      Q[i][j] = Q[i][j] / Q[i][i];
    }

    // 2a. Print the resulting form
    /*
      cout << endl;
      cout << "2a. " << endl;
      cout << "[ " << Q[1][1] << " " << Q[1][2] << " " << Q[1][3] << " ]" << endl;
      cout << "[ " << Q[2][1] << " " << Q[2][2] << " " << Q[2][3] << " ]" << endl;
      cout << "[ " << Q[3][1] << " " << Q[3][2] << " " << Q[3][3] << " ]" << endl;
      cout << endl;
    */
    
    // 3. Main Loop
    for (long k=i+1; k<=n; k++)     
      for (long l=k; l<=n; l++)     
	Q[k][l] = Q[k][l] - Q[k][i] * Q[i][l];     
    
    
    // 3a. Print the resulting form
    /*
      cout << endl;
      cout << "3a. " << endl;
      cout << "[ " << Q[1][1] << " " << Q[1][2] << " " << Q[1][3] << " ]" << endl;
      cout << "[ " << Q[2][1] << " " << Q[2][2] << " " << Q[2][3] << " ]" << endl;
      cout << "[ " << Q[3][1] << " " << Q[3][2] << " " << Q[3][3] << " ]" << endl;
      cout << endl;
    */
    
    i = i + 1;
  }


  // 4. Copy the answer to the "Cholesky" array
  counter = 0;
  for (long i=1; i<=n; i++)
    for (long j=i; j<=n; j++) {
      Cholesky[counter] = Q[i][j]; 
      counter++;
    }


  // 4a. Print the resulting form
  /*
  cout << endl;
  cout << "4a.  The Cholesky Decomposition is: " << endl;
  cout << "[ " << Q[1][1] << " " << Q[1][2] << " " << Q[1][3] << " ]" << endl;
  cout << "[ " << Q[2][1] << " " << Q[2][2] << " " << Q[2][3] << " ]" << endl;
  cout << "[ " << Q[3][1] << " " << Q[3][2] << " " << Q[3][3] << " ]" << endl;
  cout << endl;
  */


}








/////////////////////////////////////////////////////////////
// This is a modified version which takes a Matrix_mpz. =) //
// (Which is assumed to be global, hence twice the form it refers to.)
/////////////////////////////////////////////////////////////

void CholeskyDecomposition(Matrix_mpz QQ, double Cholesky[]) {

  // Check the matrix is square, and determine its dimension.
  if ((QQ.NumRows() != QQ.NumCols()) || (QQ.NumRows() <= 0)) {
    cout << "Error in CholeskyDecomposition():  The matrix is empty or not square!" << endl;
    exit(1);
  }


  long n = QQ.NumRows();
  double Q[n+1][n+1];  // Only the positive indices (1--3) are used... =)

  // 1. Initialize (from a symmetric matrix QQ) -- using 2 * QQ right now...
  long counter = 0;
  for (long i=1; i<=n; i++)
    for (long j=i; j<=n; j++) {
      Q[i][j] = 0.5 * QQ(i,j).get_d();
      counter++;
    }
  long i = 0;


  // 1a. Print the resulting form
  /*
  cout << endl;
  for (long r=1; r<=n; r++) {
    cout << "[ ";
    for (long s=1; s<=n; s++) 
      cout << Q[r][s] << " ";
    cout << "]" << endl;
  }
  cout << endl;
  */
  
  
  // 2. Loop on i
  i = i + 1;
  while (i != n) {
    for (long j=i+1; j<=n; j++) {    
      Q[j][i] = Q[i][j];
      Q[i][j] = Q[i][j] / Q[i][i];
    }

    // 2a. Print the resulting form
    /*
    cout << endl;
    for (long r=1; r<=n; r++) {
      cout << "[ ";
      for (long s=1; s<=n; s++) 
	cout << Q[r][s] << " ";
      cout << "]" << endl;
    }
    cout << endl;
    */
    
    // 3. Main Loop
    for (long k=i+1; k<=n; k++)     
      for (long l=k; l<=n; l++)     
	Q[k][l] = Q[k][l] - Q[k][i] * Q[i][l];     
    
    
    // 3a. Print the resulting form
    /*
    cout << endl;
    for (long r=1; r<=n; r++) {
      cout << "[ ";
      for (long s=1; s<=n; s++) 
      cout << Q[r][s] << " ";
      cout << "]" << endl;
    }
    cout << endl;
    */
    
    i = i + 1;
  }


  // 4. Copy the answer to the "Cholesky" array
  counter = 0;
  for (long i=1; i<=n; i++)
    for (long j=i; j<=n; j++) {
      Cholesky[counter] = Q[i][j]; 
      counter++;
    }


  // 4a. Print the resulting form
  /*
  cout << endl;
  cout << " We obtain the Cholesky Decomposition of:" << endl;
  for (long r=1; r<=n; r++) {
    cout << "[ ";
    for (long s=1; s<=n; s++) 
      cout << Q[r][s] << " ";
    cout << "]" << endl;
  }
  cout << endl;
  */


}

