
#if !defined(REPRESENTABILITY_H)
#define REPRESENTABILITY_H


// Some useful/necessary libraries
#include "GMP_class_extras/mpz_class_extras.h" 

#include "local_condition_class.h"

#include "local_split_coverings_class.h"    // NOTE: This replaces the single-cover verison...

#include "local_checking_stats.h"    // NOTE: This replaces the single-cover verison...

//#include "Boolean_Ternary_Theta_Class/boolean_ternary_theta_class.h"

#include "qf_project_class.h"
#include "qf_datafiles_class.h"  // I think this last one is unnecessary now, since it's included (and needs) QF_Project


// -----------------------------------------------------


class Representability {


 public:
  Representability(const Matrix_mpz & QQ, const double & cusp_const = -1.0, const string & type="minimal");   // Constructor
  Representability(const QF_Datafiles & qf_data, const Matrix_mpz & QQ, const double & cusp_const = -1.0, const string & type="minimal");   // Constructor
 ~Representability();                      // Destructor


  vector<mpz_class> GetExceptions() const;  // Returns the exceptions of Q.


  // Variables:
  // ----------
  bool overflow;




 private:

  // Variables:
  // ----------
  long DEFAULT_CHECKING_DEPTH;               // Sets the Default Checking depth we guarantee for checking numbers -- used in _ThetaPrecision()
  float APPROXIMATE_THETA_BOX_SIZE;          // Sets the approximate theta function box size used in _Get_Cover_Ternary_Thetas()
                                             // TO FIX: THIS SHOULD BE A LONG THROUGHOUT, NOT A FLOAT!  

  LocalSplitCoverings _local_cover;  // The local split covering using ternaries with small splitting values

  Matrix_mpz _BigForm;   // The (global) form we're checking the representability of
  long _level;          // The level of the big form
  long _character;      // The character of the big form  (i.e., the number c giving the real Dirichlet character chi(*) = (c/*) ).


  double _passed_cusp_const;     // The cusp constant passed in from the constructor...  (defaults to -1.0 if no argument is passed.)
  double _F4_bound;     // The F4 bound for representability of a number by Q 
  
  vector<long> _plist;          // List of all eligible primes
  vector<double> _f4list;       // and their associated F4(p) values.
  long _max_num_prime_factors;   // This is the maximal number of prime factors which can appear in an eligible squarefree number. (It's >=1).

  PowerSeries<mpz_class> _big_theta;  // The power series theta function of the 4-variable form.

  vector< vector<boolean_theta> > _ternary_theta;  // The boolean theta functions for the local split covering (ternary) forms.
                                                   // (THIS IS FROM THE OLDER CHECKING ROUTINES)  

  // 2 variables for the new checking routine
  vector<boolean_theta> _current_cover_ternary_thetas;  // This holds the boolean theta functions for the current local split cover.
  long _current_cover_num;    // Says the number (from 1) of the current local cover whose ternaries are stored in _current_cover_ternary_thetas.


  vector<mpz_class> _squarefree_exceptions;  // Lists the squarefree exceptions of the form

  vector<mpz_class> _squarefree_overflows;  // Lists the numbers which cause an overflow when checking for squarefree exceptions

  vector<mpz_class> _all_exceptions;  // Lists the squarefree exceptions of the form


  bool _used_squarefree_exception_file;  // Says whether we read our exceptions from a file or not... (thereby avoiding any computations!)


  long _BIGFORM_THETA_PRECISION;    // This is the precision of the quaternary theta function.  (Please no larger than 10,000)!

  vector<LocalCheckingStats> _local_stats_vec;  // Lists the checking stats for each local cover.



  // Methods: 
  // --------
  void _MakeWithMinimalCover(const QF_Datafiles & qf_data, const Matrix_mpz & QQ, const double & cusp_const);   // Constructor
  void _MakeWithUniformCover(const QF_Datafiles & qf_data, const Matrix_mpz & QQ, const double & cusp_const);


  // TO DO: Fix these next 4/5 for LocalSplitCovering(S)
  double _thetatime_for_minimal_local_cover(const LocalSplitCoverings & cover) const;
  double _thetatime_for_uniform_local_cover(const LocalSplitCoverings & cover, const mpz_class & E) const;
  mpz_class _uniform_maximum_exception_target(const LocalSplitCoverings & cover) const;

  string _choose_checking_method(const QF_Datafiles & qf_data, const double & F4) const;

  /*
  set<mpz_class> _Make_Uniform_Cover_Exceptions(const LocalSplitCoverings & cover, 
						const mpz_class & desired_exception_bound,
						const QF_Datafiles & qf_data) const;
  */


  void _GetF4Bound(const QF_Datafiles & data);           // TO DO


  void _GetEligiblePrimes(const QF_Datafiles & qf_data);
  bool _ReadEligiblePrimes(const string & primefilename, const string & F4filename);
  bool _WriteEligiblePrimes(const string & primefilename, const string & F4filename) const;
  double _F4PrimeIsoExact(const long p) const;
  void _ComputeEligiblePrimes();  
  void _FindMaxNumPrimeFactors();

  double _Upper_bound_for_eligible_numbers() const;
  mpz_class _ThetaPrecision(const long diag_coeff, const long iteration_depth = 1) const;  
  bool _MakeTernaryThetas(const QF_Datafiles & qf_data, const long iteration_depth = 1);
  void _UpdateTernaryThetas(const QF_Datafiles & qf_data);

  void _GetQuaternaryTheta(const QF_Datafiles & qf_data);
  bool _ReadQuaternaryTheta(const string & quaternary_theta_filename);

  void _CheckSquarefree(const QF_Datafiles & qf_data);
  void _WriteSquarefree(const QF_Datafiles & qf_data) const;
  bool _ReadSquarefree(const QF_Datafiles & qf_data);
  //  void _NewestCheckNumbers_by_size(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const int Max_Depth);
  void _NewestCheckNumbers_by_number(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const int Max_Depth);
  void _NewestCheckNumbers_Subroutine(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const int Max_Depth, 
				      vector<mpz_class> & possible_exception_list);


  // These are the new checking routines in "representability--new_checksquarefree.cc"
  void _CheckSquarefree_New(const QF_Datafiles & qf_data);
  void _NewestCheckNumbers_by_cover_Subroutine(const QF_Datafiles & qf_data, const vector<mpz_class> & num, 
					       const long & num_of_prime_factors, vector<mpz_class> & possible_exception_list,
					       vector<mpz_class> & possible_exception_depths);
  void _NewestCheckNumbers_by_cover(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const long & num_of_prime_factors);
  void _Get_Cover_Ternary_Thetas(const QF_Datafiles & qf_data, const long & cover_num);





  void _CheckSquareFactors();   // TO DO

};










// ---------------------------------------------------

/*
#include "readmagmaseries.cc"  // This is needed for now because of line 217 in scraps.cc...
#include "maketheta_small.cc"  // This is needed for the ThetaPrecision routine.
*/

#include "Utilities/string_utils.h"

/*  THIS IS NOT USED FOR NOW!
#include "Ternary_Exceptions/Ternary_Exceptions.h"
*/
#include "Quaternary_Exceptions/Quaternary_Exceptions.h"


// This is the real deal
#include "representability.cc"
#include "representability--new_checksquarefree.cc"


#endif
