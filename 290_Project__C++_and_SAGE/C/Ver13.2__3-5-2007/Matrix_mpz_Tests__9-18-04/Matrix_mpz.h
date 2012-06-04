#if !defined(MATRIX_MPZ_H)
#define MATRIX_MPZ_H

using namespace std;

#include <gmpxx.h>
#include <valarray>
#include <vector>
#include <iostream>



// My newer Matrix class! =)
// Made as a single valarray vector, whose dimensions we keep track of...

class Matrix_mpz {
private:
  valarray<mpz_class> M;
  size_t m, n;


  mpz_class & operator[](size_t ind);                  // Allow the notation M[n*i+j] 
  const mpz_class & operator[](size_t ind) const;      // Allow the notation M[n*i+j] 

  mpz_class & operator()(size_t ind);                  // Allow the notation M(n*i+j) 
  const mpz_class & operator()(size_t ind) const;      // Allow the notation M(n*i+j) 

  size_t Length() const;




public:
  // Want to initialize all numeric matrices to zero...
  Matrix_mpz(size_t r, size_t s);
  Matrix_mpz();

  // Copy Constructor -- Needed so that the valarray is resized! =)
  void operator=(const Matrix_mpz & source);

  // Comparison operator 
  bool Matrix_mpz::operator==(const Matrix_mpz & source) const;


  size_t NumRows() const;
  size_t NumCols() const;


  // Allow the notation M % R 
  Matrix_mpz operator%(const mpz_class & R) const;

  // Quickly evaluate the expression v^t * M * v (mod R)
  mpz_class EvaluateQuadratic(const valarray<mpz_class> & v, const mpz_class & R) const;


  mpz_class & operator()(size_t row, size_t col);                // Allow the notation M(i,j) 
  const mpz_class & operator()(size_t row, size_t col) const;    // Allow the notation M(i,j) 

  /*
  // Can't do this since it takes only one argument!
  mpz_class & operator[](size_t row, size_t col);                // Allow the notation M[i,j] 
  const mpz_class & operator[](size_t row, size_t col) const;    // Allow the notation M[i,j] 
  */


  // Define Matrix multiplication
  Matrix_mpz operator*(const Matrix_mpz & B);



  void Print(ostream & out) const;

  void PrintM(ostream & out) const;


  // Check this: would like to zero it out (now)
  //   or to preserve the existing matrix... (later?)
  void SetDims(size_t r, size_t s);

  bool IsSquare() const;

  bool IsSymmetric() const;



  // Converts a Vector_mpz into a diagonal Matrix_mpz
  void DiagonalMatrix(const valarray<mpz_class> & v);

  void Transpose();

  Matrix_mpz GetTranspose() const;

  // Evaluates T^t * Q * T for a matrix T
  Matrix_mpz EvaluateQuadratic(const Matrix_mpz & T) const;

  // Evaluates T^t * Q * T for a valarray T
  mpz_class EvaluateQuadratic(const valarray<mpz_class> & T) const;

  // Evaluates T^t * Q * T for a vector T
  mpz_class EvaluateQuadratic(const vector<mpz_class> & T) const;

  // Evaluates T^t * Q * T for a vector T (of longs)
  mpz_class EvaluateQuadratic(const vector<long> & T) const;

  // Extracts a square matrix according to the entries in Index 
  // (Note: The indexing in Index starts at 1, not at zero.) 
  Matrix_mpz ExtractSquareSubmatrix(const valarray<size_t> & Index) const;

  // Extracts a column vector from a matrix 
  valarray<mpz_class> ExtractColumn(size_t j) const;


  // Extracts a row vector from a matrix 
  valarray<mpz_class> ExtractRow(size_t i) const;


  // Extracts a square submatrix with indices associated to the given vector
  Matrix_mpz ExtractSquareSubmatrixOrdered(const valarray<size_t> & Index) const;                                    







  void SwapRows(size_t i, size_t j);   // Swaps the ith and jth rows of our matrix




  mpz_class Determinant() const;       // Computes the Determinant

  Matrix_mpz Adjoint() const;          // Computes the Adjoint

  //  Matrix_mpz Inverse()  //  <---- Need an mpq_class type to do this...

  mpz_class QFLevel() const;   // WARNING: FORGOT THE FACTOR OF 2
};




// Define the << operator for the Matrix_mpz type
ostream & operator<<(ostream & out, const Matrix_mpz & matr);

#endif
