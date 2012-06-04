//////////////////////////////////////////////////////////////////////////
//
// To Do: (9/24/04)
//   - Replace the determinant routine with something working in any ring.
//   - Make this a templated class
//   - Change EvaluateQuadratic(v,m) to EvaluateQuadraticMod(v,m)
//   - Allow vector and valarray inputs in our routines
//
//////////////////////////////////////////////////////////////////////////


#include "Matrix_mpz.h"

using namespace std;

#include <gmp.h>
#include <gmpxx.h>
#include <valarray>
#include <iostream>

#include "LCM.h"


Matrix_mpz::Matrix_mpz(size_t r, size_t s) {
    M.resize(r*s);
    M = mpz_class(0);
    m = r;
    n = s;
}

Matrix_mpz::Matrix_mpz() {
  M.resize(0);
    m = 0;
    n = 0;
}


///////////////////////////////////////////////////////////////////
// Copy Constructor -- Needed so that the valarray is resized! =)
///////////////////////////////////////////////////////////////////
void Matrix_mpz::operator=(const Matrix_mpz & source) {
  // Protect against self-assignment
  if (this != &source) {
    m = source.m;
    n = source.n;
    
    /*
    // Redid this since I heard there was a valarray resize bug!
    valarray<mpz_class> N(m*n,0);
    for(size_t i=0; i < m*n; i++)
    N[i] = source[i];
    M = N;
    */
    
    //      /*
    // This was the original way...
    M.resize(m*n);
    M = source.M;
    //      */
  }
}



// Comparison operator 
bool Matrix_mpz::operator==(const Matrix_mpz & source) const {
  
  // Check the sizes are the same
  if ((m != source.m) || (n != source.n))
    return false;
  
  // Check the entries are the same
  for(long i=1; i<=m; i++)
    for(long j=1; j<=n; j++)
      if (source(i,j) != (*this)(i,j))
	return false;
  
  // If so, then return true
  return true;
}



/* UNABLE TO CORRECTLY DEFINE THE ADDITION OF TWO MATRICES!!
   
// Define the addition of two matrices
//  void operator+(const Matrix_mpz B, const Matrix_mpz C) {
//    Matrix_mpz operator+(const Matrix_mpz B, const Matrix_mpz C) {
//    Matrix_mpz operator+(const Matrix_mpz & B, const Matrix_mpz & C) {
//    Matrix_mpz & operator+ (Matrix_mpz B, Matrix_mpz C) {
//  void operator+(Matrix_mpz & A, const Matrix_mpz & B, const Matrix_mpz & C) {
//    void operator+(Matrix_mpz & B, const Matrix_mpz & C) {
//  void operator+(Matrix_mpz B, Matrix_mpz C) {

// Check sizes are the same before adding
if ((B.m == C.m) && (B.n == C.n)) {
Matrix_mpz A; 
A.m = B.m;
A.n = B.n;

A.M.resize(m*n);
A.M = B.M + C.M;
return A;
} 
else 
cout << "Error in matrix addition: They aren't the same size!" << endl;

}

*/


size_t Matrix_mpz::NumRows() const {
  return m;
}    

size_t Matrix_mpz::NumCols() const {
  return n;
}    

size_t Matrix_mpz::Length() const {
  return m*n;
}    



//////////////////////////////
// Allow the notation M % R 
///////////////////////////////
Matrix_mpz Matrix_mpz::operator%(const mpz_class & R) const {
  
  // Check the modulus is > 0
  if (R < 0) {
    cout << "Error in % operator: The modulus must be positive!" << endl;
    exit(1);
  }
  
  Matrix_mpz M_mod(m,n);
  
  for(size_t i=1; i<=m; i++)
    for(size_t j=1; j<=n; j++)
      M_mod(i,j) = (((*this)(i,j) % R) + R) % R;
  
  return M_mod;
}




////////////////////////////////
// Allow the notation M[n*i+j] 
////////////////////////////////
mpz_class & Matrix_mpz::operator[](size_t ind) {
  return M[ind];
}

const mpz_class & Matrix_mpz::operator[](size_t ind) const {
  return M[ind];
}


/////////////////////////////////
// Allow the notation M(n*i+j) 
////////////////////////////////
mpz_class & Matrix_mpz::operator()(size_t ind) {
  return M[ind - 1];
}

const mpz_class & Matrix_mpz::operator()(size_t ind) const {
  return M[ind - 1];
}




///////////////////////////////
// Allow the notation M(i,j) 
///////////////////////////////
mpz_class & Matrix_mpz::operator()(size_t row, size_t col) {
  if ((1 <= row) && (1 <= col) && (row <= m) && (col <= n))   // This does some basic error checking
    return M[n*(row-1) + (col - 1)];
  
  
  cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
  abort();
}
const mpz_class & Matrix_mpz::operator()(size_t row, size_t col) const {
  if ((1 <= row) && (1 <= col) && (row <= m) && (col <= n))   // This does some basic error checking
    return M[n*(row-1) + (col - 1)];
  
  
  cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
  abort();
}



/*
// WARNING: We can't do this since it takes only one argument!!!

// Allow the notation M[i,j] 
  mpz_class & Matrix_mpz::operator[](size_t row, size_t col) {
    if ((0 <= row < m) && (0 <= col < n))   // This does some basic error checking
      return M[n * row + col];

    cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
    abort();
  }

  // Allow the notation M[i,j] 
  const mpz_class & Matrix_mpz::operator[](size_t row, size_t col) const {
    if ((0 <= row < m) && (0 <= col < n))   // This does some basic error checking
      return M[n * row + col];

    cerr << "Error in matrix read: row " << row << " and col " << col  << " are out of range... =|" << endl;
    abort();
  }
*/



// Define Matrix multiplication
Matrix_mpz Matrix_mpz::operator*(const Matrix_mpz & B) {
  
  // Abort if we can't multiply them
  if (NumCols() != B.NumRows()) {
    cout << "Error in Matrix Multiplication: Trying to multiply a " 
	 << NumRows() << " x " << NumCols() << " matrix by a " 
	 << B.NumRows() << " x " << B.NumCols() << " matrix." ;
    exit(1);
  }
  
  // Do the multiplication 
  Matrix_mpz C((*this).NumRows(), B.NumCols());
  for (long i=1; i<=C.NumRows(); i++)
    for (long j=1; j<=C.NumCols(); j++) {
      C(i,j) = 0;
      //      cout << "The (" << i << "," << j << ") entry is " << C(i,j) << endl; 
      for (long k=1; k<=(*this).NumCols(); k++)
	C(i,j) += (*this)(i,k) * B(k,j);
      //      cout << "The (" << i << "," << j << ") entry is " << C(i,j) << endl; 
    }
  
  // Return the product
  return C;
}



//////////////////////
// Prints the Matrix 
//////////////////////
void Matrix_mpz::Print(ostream & out) const {
  
  /*
    cout << " m = " << m << "  n = " << n << endl;
    cout << " length of the valarray = " << M.size() << endl;
  */
  
  for(size_t i = 1; i <= m; i++) {
    out << " [ ";
    for(size_t j = 1; j <= n; j++) {
      out << (*this)(i,j);
      if (j <= n - 1)
	out << ", ";
    }
    out << " ]" << endl;
  }
}


void Matrix_mpz::PrintM(ostream & out) const {
  Print(out);
}



// Check this: would like to zero it out (now)
//   or to preserve the existing matrix... (later?)
void Matrix_mpz::SetDims(size_t r, size_t s){ 
  M.resize(r*s);
  M = mpz_class(0);
  m = r;
  n = s;
}


bool Matrix_mpz::IsSquare() const {
  return (m == n);
}


bool Matrix_mpz::IsSymmetric() const {
  bool flag = false;
  
  if (m == n) { 
    flag = true;
    for(size_t i=1; i<=n; i++)
      for(size_t j=1; j<=m; j++)
	if ((*this)(i,j) != (*this)(j,i)) 
	  flag = false;
  }
  
  return flag;
  }




// Converts a Vector_mpz into a diagonal Matrix_mpz
void Matrix_mpz::DiagonalMatrix(const valarray<mpz_class> & v) {
  size_t n = v.size();
  SetDims(n,n);
  
  for(size_t i=1; i<=n; i++)
    (*this)(i,i) = v[i-1];
} 


void Matrix_mpz::Transpose() {
  valarray<mpz_class> N(n*m);  
  // cout << " m = " << m << "  n = " << n << endl; 
  
  for(size_t i=1; i<=n; i++)
    for(size_t j=1; j<=m; j++)
      N[m*(i-1) + (j-1)] = (*this)(j,i);
  
  M = N;
  swap(m,n);
}


Matrix_mpz Matrix_mpz::GetTranspose() const {
  
  Matrix_mpz Trans;
  Trans = *this;
  Trans.Transpose();
  
  return Trans;
}




//////////////////////////////////////////////////////////////////////////////
// Quickly evaluate the expression v^t * M * v (mod R)   <-- Should change the name to EvaluateQuadraticMod()...
//////////////////////////////////////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const valarray<mpz_class> & v, const mpz_class & R) const {
  mpz_class total;
  total = 0;
  
  // Check the modulus is > 0
  if (R < 0) {
    cout << "Error in EvaluateQuadratic(): The modulus must be positive!" << endl;
    exit(1);
  }
  
  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }
  
  // Check the valarray has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }
  
  // Evaluate v^T * Q * v (mod R)
  for(size_t i=1; i<=m; i++)
    for(size_t j=1; j<=n; j++)
      total = (total + (v[i-1] * (*this)(i,j) * v[j-1])) % R;
  
  /*
    if (total >= R) 
    cout << " Error in EvaluateQuadratic: R exceeded for vector: " << v << endl;
    
    if (total < 0) 
    cout << " Error in EvaluateQuadratic: Negative value for vector: " << v << endl;
  */
  
  return (total + R) % R;  // This is necessary for some strange GMP reason...  (see 2/18/04 Notes...)
}



/////////////////////////////////////////
// Evaluates T^t * Q * T for a matrix T
////////////////////////////////////////
Matrix_mpz Matrix_mpz::EvaluateQuadratic(const Matrix_mpz & T) const {

  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }
  
  // Check that the new matrix T has the correct number of rows
  if (T.NumRows() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The new matrix has the wrong number of rows!" << endl;
    exit(1);
  }
  
  Matrix_mpz ans(T.NumRows(), T.NumRows());
  
  for(size_t i=1; i<=T.NumRows(); i++)
    for(size_t j=1; j<=T.NumRows(); j++)
      for(size_t k=1; k<=m; k++)
	for(size_t l=1; l<=n; l++)
	  ans(i,j) += T(k,i) * (*this)(k,l) * T(l,j); 
  
  return ans;
}
    

/////////////////////////////////////////////
// Evaluates v^t * Q * v for a valarray v
/////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const valarray<mpz_class> & v) const{
    
  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }
  
  // Check the valarray has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }
  // Should check that *this is a square matrix, and T has the same # of rows.
  
  mpz_class ans;
  ans = 0;
  
  for(size_t k=1; k<=m; k++)
    for(size_t l=1; l<=n; l++)
      ans += v[k-1] * (*this)(k,l) * v[l-1]; 
  
  return ans;
}


///////////////////////////////////////////////////////
// Evaluates v^t * Q * v for a vector v (of mpz_class)
///////////////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const vector<mpz_class> & v) const{
   
  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }
  
  // Check the vector has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }
  
  mpz_class ans;
  ans = 0;
  
  for(size_t k=1; k<=(*this).NumRows(); k++)
    for(size_t l=1; l<=(*this).NumCols(); l++)
      ans += v[k-1] * (*this)(k,l) * v[l-1]; 
  
  return ans;
}



////////////////////////////////////////////////////
// Evaluates v^t * Q * v for a vector v (of longs) 
/////////////////////////////////////////////////////
mpz_class Matrix_mpz::EvaluateQuadratic(const vector<long> & v) const{
  
  // Check the matrix is square
  if ((*this).NumRows() != (*this).NumCols()) {
    cout << "Error in EvaluateQuadratic(): The matrix is not square!" << endl;
    exit(1);
  }
  
  // Check the vector has the appropriate size
  if (v.size() != (*this).NumRows()) {
    cout << "Error in EvaluateQuadratic(): The valarray is the wrong size!" << endl;
    exit(1);
  }
  
  mpz_class ans;
  ans = 0;
  
  for(size_t k=1; k<=(*this).NumRows(); k++)
    for(size_t l=1; l<=(*this).NumCols(); l++)
      ans += v[k-1] * (*this)(k,l) * v[l-1]; 
  
  return ans;
}





// Extracts a square matrix according to the entries in Index 
// (Note: The indexing in Index starts at 1, not at zero.) 
Matrix_mpz Matrix_mpz::ExtractSquareSubmatrix(const valarray<size_t> & Index) const {
  size_t len;
  len = Index.size(); 
  
  size_t i, j, max_ind; 
  
  // Check that the biggest entry of Index is in range of the matrix M
  max_ind = Index.max();  
  if ( (max_ind > NumRows()) || (max_ind > NumCols()) ) {
    cout << "\n Error in ExtractSquareSubmatrix: The index vector is out of range! \n\n";
    return Matrix_mpz();
  }
  
  // Extract the appropriate entries into Mnew
  Matrix_mpz Mnew(len, len);
  
  /*
    cout << " Just created Mnew: " << endl;
    Mnew.PrintM(cout);
    cout << endl;
  */
  
  for (i=1; i<=len; i++) 
    for (j=1; j<=len; j++) 
      Mnew(i,j) = (*this)(Index[i-1],Index[j-1]);
  
  /*
    cout << " Now we have Mnew: " << endl;
    Mnew.PrintM(cout);
    cout << endl;
  */
    
  return Mnew;
}



// Extracts a column vector from a matrix 

valarray<mpz_class> Matrix_mpz::ExtractColumn(size_t j) const {
  
  valarray<mpz_class> col;
  size_t i;
  
  if ( (j > n) || (j <= 0) ) {
    cout << "\n Error in ExtractColumn: The column index " \
	 << j << " exceeds the number of columns " << m << ".\n";
    return col;
  }
  
  col.resize(n);
  for(i=1; i<=m; i++)
    col[i-1] = (*this)(i,j);
  
  return col;
}



// Extracts a row vector from a matrix 

valarray<mpz_class> Matrix_mpz::ExtractRow(size_t i) const {
  
  valarray<mpz_class> row;
  size_t j;
  
  if ( (i > m) || (i <= 0) ) {
    cout << "\n Error in ExtractColumn: The row index " \
	 << i << " exceeds the number of rows " << n << ".\n";
    return row;
  }
  
  row.resize(m);
  for(j=1; j<=n; j++)
    row[j-1] = (*this)(i,j);
  
  return row;
}



// Extracts a square submatrix with indices associated to the given vector

Matrix_mpz Matrix_mpz::ExtractSquareSubmatrixOrdered(const valarray<size_t> & Index) const
{
    size_t len;
    len = Index.size();
    
    size_t i, j, max;
    // Check that the biggest entry of Index is in range of the matrix M
    max = 0;
    for (i=1; i<=len; i++)
      if (Index[i-1] > max)
	max = Index[i-i];
    
    /*
      cout << "\n ExtractSquareSubmatrixOrdered Status: \n";
      cout << " Matrix M = " << M << "\n";
      cout << " Index vector = " << Index << "\n";
      cout << " maximum index vector entry= " << max << "\n";
    */
    
    if ( (max > m) || (max > n) ) {
      cout << "\n Error in ExtractSquareSubmatrixOrdered: The index vector is out of range! \n\n";
      return *this;
    }
    
    // Extract the appropriate entries
    Matrix_mpz Mnew(len,len);
    for (i=1; i<=len; i++)
      for (j=1; j<=len; j++)
	Mnew(i,j) = (*this)(Index[i-1],Index[j-1]);
        
    return Mnew;
  }
                                                                                                                     





  // Swaps the ith and jth rows of our matrix

  void Matrix_mpz::SwapRows(size_t i, size_t j){
    
    mpz_class temp;
    size_t k;

    for (k=1; k<=n; k++) {
      temp = (*this)(i,k);
      (*this)(i,k) = (*this)(j,k);
      (*this)(j,k) = temp;
    }

  }






  // Computes the Determinant
  mpz_class Matrix_mpz::Determinant() const {
	  
    // Check that the matrix is square
    if (IsSquare() && (m>0)) {
      size_t i, j, k;
      Matrix_mpz Temp;
      Temp = *this;
      mpz_class extra = 1;

      // Check for the 1 x 1 case first 
      if (m == 1)
	return Temp(1,1);

      mpz_class A_ij, A_jj;
      // Go through the columns in order
      for (j=1; j<=n-1; j++) {  
	for (i=j+1; i<=m; i++) {  
	  
	  /*
	    cout << "Before: " << endl;
	    cout << " i = " << i << "  j = " << j << endl; 
	    cout << " Extra factor = " << extra << endl; 
	    Temp.Print(cout);
	    cout << endl; 
	  */

	  // Check to see if our pivot entry is zero
	  if (Temp(j,j) == 0) {
	    for(k=j+1; k<=n; k++)
	      if (Temp(k,j) != 0) {
		Temp.SwapRows(j,k);
		extra *= -1;
		break;
	      }
	  }

	  A_ij = Temp(i,j);
	  A_jj = Temp(j,j);


	  // If a non-zero pivot is found, use it to do an elementary row operation 
	  if (Temp(j,j) != 0) {
	  
	    // Perform the linear combimation [Row i -> - Temp(i,j) * Row i + Temp(j,j) * Row j]
	    for(k=1; k<=n; k++) {
	      //    cout << " - " << Temp(i,j) << " * " <<  Temp(j,k) << " + " << Temp(j,j) << " * " << Temp(i,k) << endl;
	      
	      Temp(i,k) = - A_ij * Temp(j,k) + A_jj * Temp(i,k);
	    }
	    
	    extra *= Temp(j,j);
	  }
	    
	   /*
	    cout << "After: " << endl;
	    cout << " i = " << i << "  j = " << j << endl; 
	    cout << " Extra factor = " << extra << endl; 
	    Temp.Print(cout);
	    cout << endl; 
	   */

	}
      }
      
      // Now multiply by the determinant of the upper triangular matrix Temp

      mpz_class diag_det;
      diag_det = 1;
      for(k=1; k<=m; k++) 
	diag_det *= Temp(k,k);

      //      cout << "the det is " << (diag_det / extra) << endl;
      
      return (diag_det / extra);
    }


    cerr << " Error in Determinant method: The matrix is not square or it has size 0! =(" << endl;
    cerr << " Using Q = ";
    PrintM(cout);
    cerr << endl;
    cerr << " where m = " << m << "  and  n = " << n << endl;
    abort();
  } 


  Matrix_mpz Matrix_mpz::Adjoint() const {

    // Note: The Adjoint matrix should be square...
    Matrix_mpz Adj(m,n);

    for(size_t i=1; i<=m; i++)
      for(size_t j=1; j<=n; j++){

	size_t k_ptr, l_ptr;
	Matrix_mpz Minor_mat(m-1, n-1);

	//	cout << "\n i = " << i << "  and  j = " << j << endl << endl;

	// Loop to construct the (m-1) x (n-1) submatrix
	k_ptr = 1; 
	for(size_t k=1; k<=m-1; k++) {

	  if (k == i)   // Skip the i-th row  
	    k_ptr++;

	  l_ptr = 1;	    
	  for(size_t l=1; l<=n-1; l++){

       	    if (l == j)   // Skip the j-th column
		l_ptr++;

	    /*  
	    cout << " k = " << k << "  and  l = " << l << endl;
	    cout << " k_ptr = " << k_ptr << "  and  l_ptr = " << l_ptr << endl;
	    */
	    
	    Minor_mat(k, l) = (*this)(k_ptr, l_ptr);

	    l_ptr++;
	  }

	  k_ptr++;
	}

	/*
	cout << " The matrix minor is \n"; 
	Minor_mat.PrintM(cout); 
	cout << endl;
	*/	

	mpz_class sign;
	if ((i+j) % 2 == 0) 
	  sign = mpz_class(-1);
	else 
	  sign = mpz_class(1);


	Adj(j,i) = sign * (Minor_mat).Determinant();

      }

    return Adj;

  }



  //  Matrix_mpz Inverse()  //  <---- Need an mpq_class type to do this...




  mpz_class Matrix_mpz::QFLevel() const {  // WARNING: FORGOT THE FACTOR OF 2

    mpz_class temp_lvl, det;
    det = (*this).Determinant();
    temp_lvl = 1;
    Matrix_mpz Adj;
    Adj = (*this).Adjoint();
    mpq_class temp_inv_entry;

    // Find the LCM of the denominators of the upper triangular entries
    for(size_t i=1; i<=m-1; i++)
      for(size_t j=i+1; j<=n; j++) {
	temp_inv_entry = mpq_class(Adj(i,j), det);
	temp_inv_entry.canonicalize();
	temp_inv_entry = temp_inv_entry.get_den();
	temp_lvl = LCM(temp_lvl, temp_inv_entry);
      }    

    // Find the LCM of the denominators of the diagonal entries 
    for(size_t i=1; i<=m; i++) {
      temp_inv_entry = mpq_class(Adj(i,i), det);
      temp_inv_entry.canonicalize();
      temp_inv_entry = temp_inv_entry.get_den();
      temp_lvl = LCM(temp_lvl, 4 * temp_inv_entry);
    }

    return temp_lvl;
	  
  }



// Define the << operator for the Matrix_mpz type
ostream & operator<<(ostream & out, const Matrix_mpz & matr) {    // Why doesn't "(...., Matrix_mpz & matr)" work?  =|
  matr.Print(out);
  return out;
}


