#if !defined(MATRIX_MPZ_H)
#define MATRIX_MPZ_H


using namespace std;

#include <gmpxx.h>
#include <valarray>
#include <vector>
#include <set>
#include <iostream>

#include <assert.h>

#include "../power_series.h"   // This needs to by put above the class declaration...
                               // Note: TO DO: This can replace the readseries stuff included at the end.

#include "../Utilities/string_utils.h"  
#include "../Utilities/file_utilities.h"  



// My newer Matrix class! =)
// Made as a single valarray vector, whose dimensions we keep track of...

class Matrix_mpz {

public:
  Matrix_mpz(size_t r, size_t s);   // Note: All new matrices are initalized to zero! =)
  Matrix_mpz();


  // Copy Constructor -- Needed so that the valarray is resized! =)
  void operator=(const Matrix_mpz & source);

  // Comparison operator 
  bool operator==(const Matrix_mpz & source) const;


  size_t NumRows() const;
  size_t NumCols() const;




  // Some Modulus Operators:
  // -----------------------

  // Allow the notation M % R 
  Matrix_mpz operator%(const mpz_class & R) const;

  // Quickly evaluate the expression v^t * M * v (mod R)
  mpz_class EvaluateQuadratic(const valarray<mpz_class> & v, const mpz_class & R) const;  // TO DO: We should rename this!!!






  mpz_class & operator()(size_t row, size_t col);                // Allow the notation M(i,j) 
  const mpz_class & operator()(size_t row, size_t col) const;    // Allow the notation M(i,j) 

  /*
  // Can't do this since it takes only one argument!
  mpz_class & operator[](size_t row, size_t col);                // Allow the notation M[i,j] 
  const mpz_class & operator[](size_t row, size_t col) const;    // Allow the notation M[i,j] 
  */



  // Define Matrix multiplication
  Matrix_mpz operator*(const Matrix_mpz & B);




  // Some Print Routines:
  // --------------------
  void Print(ostream & out) const;
  void PrintM(ostream & out) const;   // TO DO: Eliminate this since it's redundant!



  // Rigid I/O routines for making persistent datafiles:
  // ---------------------------------------------------
  void _FileInput(istream & in);
  void _FileOutput(ostream & out) const;   
  


  // QF Filename Generation:
  // -----------------------
  string QF_String() const;


  // Check this: would like to zero it out (now)
  //   or to preserve the existing matrix... (later?)
  void SetDims(size_t r, size_t s);    // To Do: We'd like to depricate this routine, and replace it by a SafeResize(m,n) routine...
  void SafeResize(size_t r, size_t s);    
  


  // Some Useful Tests:
  // ------------------
  bool IsSquare() const;
  bool IsSymmetric() const;
  bool IsQuadraticForm() const;


  // Converts a Vector_mpz into a diagonal Matrix_mpz
  void DiagonalMatrix(const valarray<mpz_class> & v);
  void IdentityMatrix(size_t n);  


  void Transpose();

  Matrix_mpz GetTranspose() const;




  // Evaluates Q[x] for many types: 
  // -------------------------------

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



  // Determinant, Adjoint, and Level:
  // --------------------------------
  mpz_class Determinant() const;       // Computes the Determinant
  Matrix_mpz Adjoint() const;          // Computes the Adjoint
  //  Matrix_mpz Inverse()  //  <---- Need an mpq_class type to do this...
  mpz_class QFLevel() const;   // WARNING: FORGOT THE FACTOR OF 2


  // Cholesky Decomposition:
  // -----------------------
  vector< vector<double> > CholeskyDecomposition() const;    // Gives it as a lower-triangular "matrix" and uses indices 1 -- n  (*** CHANGE THIS LATER! ***)


  // Compute the Theta function:                  <--- ADD EISENSTEIN SERIES HERE TOO!
  // ---------------------------
  PowerSeries<mpz_class> ComputeTheta(const unsigned long & precision) const;


  // ================================== Reduction.cc ================================================
  void ReduceOffDiagonal();



  // ================================== Local_Normal.cc ==================================================

  // Symmetric Elementary (Row/Col) Operations:
  // ------------------------------------------
  void SwapSymmetric(size_t i, size_t j);
  void MultiplySymmetric(mpz_class c, size_t i);
  void DivideSymmetric(mpz_class c, size_t i);
  void AddSymmetric(mpz_class c, size_t i, size_t j);

  // Local Normal form of a matrix:
  // ------------------------------
  Matrix_mpz GetLocalNormal(const mpz_class & p) const;
  


  // =================================== Local_Density_Front.cc =========================================

  // Local Density Routines:
  // -----------------------
  mpq_class Local_Density(const mpz_class & p, const mpz_class & m) const;
  mpq_class Local_Primitive_Density(const mpz_class & p, const mpz_class & m) const;


  // =================================== Local_Density_Front.cc =========================================

  // Local Invariant Routines:
  // -------------------------
  Matrix_mpz LocalDiagonal(const mpz_class & p) const;
  long HasseInvariant(const mpz_class & p) const;
  bool IsAnisotropic(const mpz_class & p) const;
  bool IsIsotropic(const mpz_class & p) const;
  valarray<mpz_class> AnisotropicPrimes() const;

  Matrix_mpz RationalDiagonal() const;   // This one may never be used... so we should probably can it for now...


  // =================================== Local_Constants.cc =========================================

  // Local Constants Routines:
  // -------------------------
  mpq_class LocalConstantCp(const mpz_class & p, const mpz_class & T) const;
  bool IsStable(const mpz_class & p, const mpz_class & T) const;


  // =================================== Siegel_Product.cc ==========================================
  mpq_class SiegelProduct(mpz_class u) const;
  void CheckSiegelRange(const long first, const long last, const string & eisfilename = "") const;

  void Check_ComputeTheta_vs_MagmaTheta() const;                                               // <------------- This is mis-filed here...

  // =================================== Local_Diagnostic.cc ========================================
  bool CheckSiegel(const mpz_class & m, const PowerSeries<mpq_class> & Eis) const;
  void CheckLocalDensity(const mpz_class & p, const mpz_class & m) const;


  // =================================== Eisenstein_Bound.cc ========================================
  mpq_class GetEisensteinLowerBound(const string & Eis_Bound_Dir) const;
  mpq_class NumericalEisensteinLowerBound(const string & Eis_Dir, const unsigned long precision = 1000) const;

  PowerSeries<mpq_class> GetMagmaEisSeries(const string & Eis_Dir, const unsigned long precision) const;
  PowerSeries<mpz_class> GetMagmaThetaSeries(const string & Theta_Dir, const unsigned long precision) const;



private:
  vector<mpz_class> M;
  size_t m, n;

  // Some silly, implementation-dependent things:
  // --------------------------------------------
  mpz_class & operator[](size_t ind);                  // Allow the notation M[n*i+j] 
  const mpz_class & operator[](size_t ind) const;      // Allow the notation M[n*i+j] 

  mpz_class & operator()(size_t ind);                  // Allow the notation M(n*i+j) 
  const mpz_class & operator()(size_t ind) const;      // Allow the notation M(n*i+j) 

  size_t Length() const;



  // =================================== Local_Density_Front.cc =========================================

  // Internal Local Density Routines:
  // --------------------------------
  mpq_class Local_Good_Density(const mpz_class & p, const mpz_class & m) const;
  mpq_class Local_Zero_Density(const mpz_class & p, const mpz_class & m) const;
  mpq_class Local_Bad_Density(const mpz_class & p, const mpz_class & m) const;
  mpq_class Local_BadI_Density(const mpz_class & p, const mpz_class & m) const;
  mpq_class Local_BadII_Density(const mpz_class & p, const mpz_class & m) const;

  
  // =================================== Local_Density_Congruence.cc =========================================

  valarray<size_t> ReindexVectorFromExtraction(const valarray<size_t> & Original, const valarray<size_t> & Extracted) const;

  mpz_class GaussLocal(size_t n, const mpz_class & p, const mpz_class & m, const mpz_class & Qdet) const;

  mpq_class Local_Good_Density_Congruence_Odd(const mpz_class & p, const mpz_class & m,
					      const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;
  mpq_class Local_Good_Density_Congruence_Even(const mpz_class & p, const mpz_class & m,
					       const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;

  mpq_class Local_Good_Density_Congruence(const mpz_class & p, const mpz_class & m,
					  const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;
  mpq_class Local_Zero_Density_Congruence(const mpz_class & p, const mpz_class & m,
					  const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;
  mpq_class Local_BadI_Density_Congruence(const mpz_class & p, const mpz_class & m,
					  const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const; 
  mpq_class Local_BadII_Density_Congruence(const mpz_class & p, const mpz_class & m,
					   const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;
  mpq_class Local_Bad_Density_Congruence(const mpz_class & p, const mpz_class & m,
					 const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;

  mpq_class Local_Density_Congruence(const mpz_class & p, const mpz_class & m, 
				     const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;
  mpq_class Local_Primitive_Density_Congruence(const mpz_class & p, const mpz_class & m,
					       const valarray<size_t> & Zvec, const valarray<size_t> & NZvec) const;



  // ========================================= CountLocal.cc ==================================================


  void Increment(valarray<mpz_class> & v, const mpz_class & R) const;

  mpz_class CountLocalNaive(const mpz_class & m, const mpz_class & R) const;

  valarray <mpz_class> CountLocalNaiveValues(const mpz_class & R) const;

  bool IsLocalSolutionType(const mpz_class & p, const valarray<mpz_class> & w, size_t solntype, 
			   const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;

  mpz_class CountLocalTypeNaive(const mpz_class & p, unsigned long k, const mpz_class & m, size_t solntype, 
				const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;



  // ========================================= CountLocal2.cc =================================================


  size_t GetLocalSolutionType(const mpz_class & p, const valarray<mpz_class> & w, 
			      const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;

  valarray<mpz_class> CountAllLocalTypesNaive(const mpz_class & p, unsigned long k, const mpz_class & m, 
					       const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;


  mpz_class CountLocalType(const mpz_class & p, long k, const mpz_class & m, size_t solntype, 
			   const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;

  mpz_class CountLocalGoodType(const mpz_class & p, long k, const mpz_class & m, 
			   const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;

  mpz_class CountLocalZeroType(const mpz_class & p, long k, const mpz_class & m, 
			   const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;
  
  mpz_class CountLocalBadType(const mpz_class & p, long k, const mpz_class & m, 
			   const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;
  
  mpz_class CountLocalBadTypeI(const mpz_class & p, long k, const mpz_class & m, 
			   const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;
  
  mpz_class CountLocalBadTypeII(const mpz_class & p, long k, const mpz_class & m, 
			   const valarray<size_t> & zero, const valarray<size_t> & nonzero) const;



  // =================================== Eisenstein_Bound.cc ========================================


  valarray<mpz_class> _PrimSqRepns(const mpz_class & p) const;  

  valarray<mpq_class> _C4_squareclass_constants(const Matrix_mpz & Q, const mpz_class & p, const valarray<mpz_class> & sqclasslist) const;

  mpq_class _adjustment_secret(const mpz_class & m, const mpz_class & N, const mpz_class & chi_top) const;

  mpz_class _integer_part_away_from_set(const mpz_class & m, const valarray<mpz_class> & S) const;

  mpq_class _Lambda4_test(const Matrix_mpz & Q, const PowerSeries<mpq_class> & EE) const;

  mpq_class _new_C4_rational(const Matrix_mpz & Q) const;



  // =============================================================================================================

  // These are deprecated... =(
  PowerSeries<mpq_class> _GetEisData(const string & Eis_Dir, const unsigned long precision) const;
  void _ComputeEisData(const string & EisFilename, const unsigned long Eis_precision) const;

  // Newer routine... =)
  string _GetMagmaComputation(const string & Filename_fragment, const string & Client_Scriptname, const string & var_string) const;


};




// Define the << operator for the Matrix_mpz type
ostream & operator<<(ostream & out, const Matrix_mpz & matr);


#include "../LCM.h"   // Do we need this??
#include "../Utilities/string_utils.h"    // Needed for the filenames in Eisenstein_Bound.cc

#include "../output.h"   // Do we need this??


// Basic Matrix Stuff
#include "Matrix_mpz.cc"


// Some Extras we need
#include "../GMP_class_extras/mpz_class_extras.h"
#include "../GMP_class_extras/vectors.h"
#include "../GMP_class_extras/mpq_stuff.h"
#include "../GMP_class_extras/Bernoulli.h"

#include "../Siegel_Diagnostic/siegel_diagnostic.h"   // Needed for the ReadSeries() routine used in "Eisenstein_Bound.cc"

// Reduction
#include "Reduction.cc"


// Local Stuff 
#include "Local_Normal.cc"
#include "Matrix_mpz_Extras.cc"   

#include "Local_Density_Front.cc" 
#include "CountLocal.cc"
#include "CountLocal2.cc"
#include "Local_Density_Congruence.cc" 

#include "Local_Invariants.cc"
#include "Local_Constants.cc"

#include "Siegel_Product.cc"
#include "Local_Diagnostic.cc"

#include "Eisenstein_Bound.cc"


#endif 
