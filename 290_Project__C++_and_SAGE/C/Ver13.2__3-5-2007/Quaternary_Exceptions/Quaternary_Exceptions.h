#if !defined(QUATERNARY_EXCEPTIONS_H)
#define QUATERNARY_EXCEPTIONS_H



class QuaternaryExceptions {
  
 public:
  QuaternaryExceptions(const Matrix_mpz & T, const mpz_class & precision,   // Computes the exceptions up to a given precision.
		    const string & quaternary_exceptions_dir = "");   
  
  QuaternaryExceptions(const Matrix_mpz & T, const set<mpz_class> & possible_set);   // Computes the exceptions from a given set of possibilities

  
  mpz_class GetPrecision() const;           // Give the precision x the exceptions are known up to (i.e. all numbers <= x).
  set<mpz_class> GetExceptions() const;     // Give the set of exceptions (up to _precision).

  bool VerifyExceptions(const mpz_class & Max=0) const;          // Verify the exceptions (up to the smaller of max and _precision).

  

 private:
  Matrix_mpz _T;                                // The matrix of the quaternary form 2*T.
  LocalConditions _local_conditions;            // The local conditions for T.
  set<mpz_class> _exception_set;                // Set of exceptions of the form T.
  mpz_class _precision;                         // The desired precision of the set _exception_set.

  set<mpz_class> _starting_set;
  set<mpz_class> _exceptional_subset;



  void _ReadExceptions(const string & filename,                                                  // Reads a ternary exception file
			       set<mpz_class> & new_exception_set, mpz_class & new_precision) const;
  void _WriteExceptions(const string & filename,                                                 // Writes a ternary exception file
			       const set<mpz_class> & new_exception_set, const mpz_class & new_precision) const;




  //  _Make(const Matrix_mpz & T, const mpz_class & precision);  // The real initialization routine
  mpz_class _MakeNiceLocalModulus() const;                    // Helper routine for _QuaternaryExceptionsUpTo().


  // Computes the set of exceptions "_exception_set".
  void _QuaternaryExceptionsUpTo(const mpz_class & desired_precision, const string & ternary_exceptions_dir);   

  // Checks if a given number is represented by a quaternary form
  bool _IsQuaternaryException(const mpz_class & num);

  mpz_class _Jagy_IntSqrt(const mpz_class & n) const;    // Jagy's IntSqrt routine used for his ternary exception finder. =)

};



// ------------------------------

#include "../output.h"               // Needed for operator<< for set<T> in VerifyExceptions(). =)
                                     // Also Needed for ReadSet_mpz_class in Read_Quaternary_Exceptions().


//#include "../Utilities/set_utils.h"  // Also Needed for ReadSet_mpz_class in Read_Quaternary_Exceptions().

#include "Quaternary_Exceptions.cc" 


#endif 
