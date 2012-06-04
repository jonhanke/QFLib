
#include<sstream>

//////////////////////////////////////////////
/// Converts a Matrix_mpz to a PARI matrix
//////////////////////////////////////////////

GEN Matrix_mpz__to__PARI(const Matrix_mpz & M) {

  //  cout << " Starting Matrix_mpz --> PARI Conversion" << endl;

  // Make the associated string stream
  ostringstream mat;
  mat << "[";
  for(long i=1; i<=M.NumRows(); i++)
    for(long j=1; j<=M.NumCols(); j++) {
      mat << M(i,j);
      /*
      cout << " (i,j) = " << i << "," << j << endl;
      cout << " M(i,j) = " << M(i,j) << endl; 
      */
      if ((i <= M.NumRows()) && (j < M.NumCols()))
	mat << ", ";
      else if ((i < M.NumRows()) && (j == M.NumCols()))
	mat << "; ";
      else 
	mat << "]";
    }
  
  // Convert that to a PARI t_MAT type!
  string str = mat.str();
  /*
  cout << " The matrix as an input string is: " << str << endl;
  */
  char* cstr = (char*) str.c_str();
  GEN PARI_MATRIX = flisexpr(cstr);

  //  cout << " Finished Matrix_mpz --> PARI Conversion" << endl;

  return PARI_MATRIX;
}

//////////////////////////////////////////////
/// Converts a PARI matrix to a Matrix_mpz --------- WARNING: THIS IS A DUMMY ROUTINE!!!
//////////////////////////////////////////////

Matrix_mpz PARI__to__Matrix_mpz(const GEN & A) {

  //  cout << " Starting PARI --> Matrix_mpz Conversion" << endl;

  /*
  // Make the associated string stream
  ostringstream mat;
  mat << "[";
  for(long i=1; i<=M.NumRows(); i++)
    for(long j=1; j<=M.NumCols(); j++) {
      mat << M(i,j);
      if ((i <= M.NumRows()) && (j < M.NumCols()))
	mat << ", ";
      else if ((i < M.NumRows()) && (j == M.NumCols()))
	mat << "; ";
      else 
	mat << "]";
    }
  */


  // Convert that to a PARI t_MAT type!
  char* cstr;
  cstr = GENtostr(A);
  string str(cstr);

  outmat(A);
  //  cout << " str = " << str << endl;
  //  cout << " Finished PARI --> Matrix_mpz Conversion" << endl;

  Matrix_mpz dummy_matrix(1,1);

  return dummy_matrix;
}


