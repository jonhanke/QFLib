

#if !defined(LOCAL_SPLIT_COVERINGS_CLASS_H)
#define LOCAL_SPLIT_COVERINGS_CLASS_H


#include "output.h"   // This is needed for the Strict I/O routines...



class LocalSplitCoverings {

public:
  LocalSplitCoverings();  // Empty Constructor ... I don't really like this, but it's needed for some future class declarations...
  LocalSplitCoverings(const string & cover_type, const Matrix_mpz & QQ);    // Constructor which calls Make(), but doesn't save the results in a file
  LocalSplitCoverings(const string & cover_type, const Matrix_mpz & QQ, 
		     const string & LocalCoverDir); // Constructor which calls Make(), but does save the results in a file
  ~LocalSplitCoverings();  // Destructor

  void Make(const string & cover_type, const Matrix_mpz & QQ, const string & LocalCoverDir);  // Constructor which finds the local conditions
  void IncrementLocalSplitCoverings();      // Adds an additional local cover to a pre-defined local covering (and writes the covering)!

  unsigned long multiplicity() const;          // Gives the number of local covers stored
  unsigned long size(const long i) const;      // Gives the number of forms involved in the i-th covering (i>=1)


  Matrix_mpz Big_Form() const;
  LocalConditions Big_Local() const;
  Matrix_mpz Get_Form(const long cover_number, const long form_number) const;
  long Get_Value(const long cover_number, const long form_number) const;
  LocalConditions Get_Local_Conditions(const long cover_number, const long form_number) const;


  // Rigid File I/O routines:              <------  STILL TO DO...
  // ------------------------
  void _FileOutput(ostream & out) const;
  void _FileInput(istream & in);


  
private:
  Matrix_mpz _big_form;            // This is the quadratic form whose local conditions we are covering, 
  LocalConditions _big_local;      // and here are it's local conditions!
  vector< vector<Matrix_mpz> > _form_list;         // Stores the associated forms
  vector< vector<long> > _d_list;                  // Stores the associated d's
  vector< vector<LocalConditions> > _local_list;   // Stores the associated local conditions for each form
  string _cover_type;              // Stores the kind of local cover

  set< vector<long> > __used_indices;     // A set of indices of the form (d,i) we have already used


  string _filename_string;         // Stores the filename we write the cover to...  
                                   // NOTE: This is *not* saved in the file, but newly generated each time by the read/make operation!


  // This is temporary data (i.e not saved/loaded) used for adding local coverings:
  // ------------------------------------------------------------------------------
  long __Largest_size_used;       // Says the length of the largest vector previously used.
  vector< vector<Matrix_mpz> > __complementary_matrix_vec;
  vector< vector<LocalConditions> > __local_conditions_vec;



  // Finds all small vectors for QQ, ordered by the numbers they represent.
  vector< vector< vector<long> > >  _ThetaVectors(const Matrix_mpz & QQ, const unsigned long bound) const;

  // Finds the (ternary) form orthogonal to the vector vec inside of BigForm.
  Matrix_mpz _FindComplementarySublattice(const Matrix_mpz & QQ, const vector<long> & vec) const;


  // Finds a local split covering of Big_Form, perpendicular to vectors of minimal length.
  void _IncrementMinimalLocalSplitCoverings();


  // --------------------------------------------------------------




  // Checks if a vector of local coverings admits single form local covers
  // (Returns empty set <==> No form is a local cover!)
  vector< set< vector<long> > >  _Find_Single_Form_Local_Covers(const vector< vector<LocalConditions> > & local_vec, 
								const vector< set< vector<long> > > & old_cover_indices,
								set< vector<long> > & excluded_indices) const;


  /*

  // -----------------------------------------------------------------------------------------------------
  // These 3 are not used... they also have some code for multiple forms, but only locally universal ones 
  // -----------------------------------------------------------------------------------------------------

  // Checks if a vector of local coverings admits a locally universal subset
  // (Returns empty set <==> No subset is a universal cover!)
  vector< set< vector<long> > >  _Find_Locally_Universal_Subset(const vector< vector<LocalConditions> > & local_vec, 
								const vector< set< vector<long> > > & old_cover_indices,
								set< vector<long> > & excluded_indices) const;

  // Checks to see which primes p have local conditions which project to give no obstruction at p
  set<mpz_class> _Find_Eligible_Covering_Primes(const vector< vector<LocalConditions> > & local_vec, 
						const set< vector<long> > & excluded_indices) const;
  
  // Finds coverings involving these primes, if it exists.
  vector< set< vector<long> > > _Find_Coverings_Using_PrimeSet(const vector< vector<LocalConditions> > & local_vec, 
					  const set< vector<long> > & excluded_indices,
					  const set<mpz_class> & prime_set) const;
  

  */



  // ---------------------------------------------------------------
  
  // Writes a local covering from a file
  void _WriteLocalSplitCoverings() const;

  // Reads a local covering from a file
  void _ReadLocalSplitCoverings();




  // ------------------------------------ Testing Routines ----------------------------------------------

  void _Testing_ThetaVectors() const;



};




// Here are some extras:
// ---------------------
//ostream & operator<<(ostream & out, const LocalSplitCovering & local);



// -------------------------------------------------------------------------------------



// Needed to find the local conditions of a form
#include "local_condition_class.h"


#include "local_split_coverings_class.cc" 
//#include "local_split_coverings_class--multiple_forms.cc" 



#endif






