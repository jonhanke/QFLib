

#if !defined(LOCAL_CONDITION_CLASS_H)
#define LOCAL_CONDITION_CLASS_H


//#include "local_score_class.h"  // Needed for storing the local scores

#include "output.h"  // Needed for output and strict I/O routines...



class LocalConditions {

public:
  LocalConditions();  // Constructor
  LocalConditions(const Matrix_mpz & QQ, const string & type = "fast_at_2");  // Constructor
  ~LocalConditions();  // Destructor

  vector<long> Get_Anisotropic_Primes() const;
  vector<long> Get_Local_Mod_Vector() const;
  vector< vector<long> > Get_Local_Repn_Array() const;   // Note: The inner vector<long> will always have size 9.

  vector<long> List_Bad_Primes() const;   // Lists primes which have a congruence obstruction
  //  vector<long> List_Bad_SquareClasses(const long p);   // Lists bad squareclasses at p.



  void IntersectWith(const LocalConditions & B);   // Replace it with the intersection of this with B 
  void UnionWith(const LocalConditions & B);       // Replace it with the union of this with B   -- THIS SHOULD NOT BE USED
  void PrimewiseUnionWith(const LocalConditions & B);  // THIS IS ONLY USED IN LocalSplitCovering::_Find_Eligible_Covering_Primes(), 
                                                       // since it's only appropriate there!

  void Make(const Matrix_mpz & QQ);  // Initializes the local conditions for a (global) quadratic form
  void Make2(const Matrix_mpz & QQ);  // Initializes the local conditions for a (global) quadratic form (this separates out the local computations)

  bool operator==(const LocalConditions & B) const;  // Compares when two local_repn_arrays are the same.
  bool operator!=(const LocalConditions & B) const;  // Compares when two local_repn_arrays are the same.

  ostream & Print(ostream & out) const;   // Prints the local conditions



  // Rigid I/O routines for making persistent datafiles (called by front-end functions for consistency):
  // ---------------------------------------------------------------------------------------------------
  void _FileOutput(ostream & out) const;
  void _FileInput(istream & in);
  

  /* NOT USED NOW...

  // Scoring Routines:
  // -----------------
  LocalScore MissedBy(const LocalConditions & B) const;  // Gives a score of local conditions of (*this) missed by B. 
  LocalScore AddedBy(const LocalConditions & B) const;   // Gives a score of local conditions of (*this) added by B. 
  */
  

  // Misc Tests:
  // -----------
  bool IsEmpty() const;         // Says if the local conditions are empty
  bool IsUniversal() const;     // Says if the local conditions are universal
  bool IsUniversalAtPrime(const mpz_class & p) const;    // Says if the local conditions are universal at the prime p
  bool IsLocallyRepresented__Squarefree(const mpz_class & m) const;  // Says whether the number is locally represented.
                                                         // Requires m to be non-negative... for now.
  bool IsLocallyRepresented(const mpz_class & m) const;  // Says whether the number is locally represented.
                                                         // Requires m to be non-negative... for now.

  
private:
  vector<long> aniso_vector;       // Gives a vector of the anisotropic primes. (Note: Destroyed under union and intersection)
  vector<long> local_mod_vector;   // Gives a vector of the highest prime-power moduli to check. (Note: Destroyed under union and intersection)
  vector< vector<long> > local_repn_array;         // Note: The inner vector<long> will always have size 9.

  bool _empty_flag;                 // This flag tells us if the local_condition is empty... =)
  bool _aniso__mod_meaning_flag;    // This flag tells us if the aniso_vector and local_mod_vector have meaning.

  /* NOT USED NOW... 
  LocalScore MissedBy_Score(const vector<long> & a, const vector<long> & b) const;
  LocalScore AddedBy_Score(const vector<long> & a, const vector<long> & b) const;
  */

  vector<long> Make_Local_Repn_Vector_At_Prime(const Matrix_mpz & local_normal_at_p, const mpz_class & NN, const mpz_class & p) const;

  void Cleanup();  // Erases irrelevant entries in the local_repn_array


  // ----------------------------------------------

  // Routines to find local conditions quickly:
  // ------------------------------------------
  vector<long> _Multiply_Local_Conditions_by_unit_at_p(const vector<long> & conditions, const mpz_class & u, const mpz_class & p) const;
  vector<long> _Multiply_Local_Conditions_by_power_at_p(const vector<long> & conditions, const long & nu, const mpz_class & p) const;
  //vector<long> _Multiply_Local_Conditions_at_2_by_power(const vector<long> & conditions, const long & nu) const;
  //vector<long> _Multiply_Local_Conditions_at_2_by_unit(const vector<long> & conditions, const mpz_class & u) const;
  vector< vector<long> > _Parse_Local_Normal_Form_at_p(const Matrix_mpz & normal_form, const mpz_class & p) const;
  vector<long> _Local_Conditions_for_elementary_form_at_p(const vector<long> & elementary_form_vector, const mpz_class & p) const;
  vector<long> _Local_Conditions_Q_at_p(const Matrix_mpz & QQ, const mpz_class & p) const;








};




// Here are some extras:
// ---------------------
ostream & operator<<(ostream & out, const LocalConditions & local);
LocalConditions Union(const vector<LocalConditions> & local_vec);
LocalConditions Intersection(const vector<LocalConditions> & local_vec);








// -------------------------------------------------------------------------------------

// TO DO:  Add a flag for representability at the real numbers: 
//         pos def = 1, neg def = -1, indef = 0.



#include "local_condition_class.cc"
//#include "local_condition_class--score_stuff.cc"   // This has all of the local score stuff


#endif






