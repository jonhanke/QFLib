

#if !defined(BOOLEAN_THETA_CLASS_H)
#define BOOLEAN_THETA_CLASS_H




#include <fstream>
#include <iostream>


/* -------------------  NOT INCLUDING THIS DISTRIBUTED STUFF FOR NOW! -----------------------

// This is the (very simple) distributed class
#include "../distributed_class.h"
*/

// This is the routine for computing Cholesky decompositions
#include "../cholesky_distributed.h"




// These are the min/max routines
#include "../max-min.h"

// This has qf_ternary_level() --- This should be moved into the qf class...
#include "../misc_utilities.h"

// This has the LCM routine --- This is also replicated in our local densities computations
#include "../LCM.h"


// This has the LCM routine --- This is also replicated in our local densities computations
//#include "Version_XXVII__9-24-04/Classes/mpz_class_extras.h"





// This is the class to deal with boolean ternary theta functions
class boolean_theta{
public:
  boolean_theta();  // Constructor
  boolean_theta(mpz_class new_precision);  // Constructor
  boolean_theta(const Matrix_mpz & Q, mpz_class new_precision);  // Constructor

  // Copy Constructor -- Needed so that the valarray is resized! =)
  void operator=(const boolean_theta & source);

  // Comparison operators 
  bool operator==(const boolean_theta & b) const;
  bool operator<=(const boolean_theta & b) const;
  bool operator>=(const boolean_theta & b) const;

  mpz_class precision() const;
  bool get_value(mpz_class i) const;
  void set_value(mpz_class i);
  void clear_value(mpz_class i);
  Matrix_mpz get_qf();
  bool read(const string & filename = "");  // Take in an optional filename...(for read and write and read_magma)
  void write(const string & filename = "") const;
  void read_magma(const string & seriesfilename);



/* -------------------  NOT INCLUDING THIS DISTRIBUTED STUFF FOR NOW! -----------------------
  void compute();
  void compute_distributed(unsigned long SLICE_MIN = 1000000);
  void compute_client(const long slice_number);

  void compute_distributed_new();
  void compute_client_new();
*/


  //  void compute_distributed_OLD(unsigned long SLICE_MIN = 1000000);

  void extended_comparison(const boolean_theta & b);



  //////////////////////////////////////////
  // Global constants used to store files //
  //////////////////////////////////////////
  
  // Says where the theta functions are
  char* THETA_DIR; 
  
  // Says where the intermediate computations are stored
  char* DIST_DIR; 

  // Says where the ternary_theta_client routine is
  //const char EXEC_DIR[] = "/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/Ver7__8-19-04/";
  char* EXEC_DIR; 
  
  // *** WARNING: This requires we setup symlinks to out favoite executables... ternary_theta_client and ternary_theta_client_new
  
  
  



  // Display Methods
  // ---------------
  // show_range(unsigned long a, unsigned long b);
  // show_first(unsigned long n);
  // show_last(unsigned long n);


  //  bool compare();



private:
  // Garrett's Diagnostic Modification
  //  boolean_ternary_theta(const boolean_ternary_theta&);



  // This should be changed to our matrix type soon...
  //
  // Also change allmethods to use _QQ instead of QQ!!!
  //
  Matrix_mpz QQ; // The quadratic form to use.  (Note: Actually, it's the matrix of 2Q!)

  vector<unsigned long> _theta;   // Stores the boolean bits of the theta function
  mpz_class _theta_precision;   // Says the highest meaningful theta coefficient  
                                    // ( We'd like this to be const, but it messes up the copy constructor!!! =( )



  unsigned long _length() const;
  void _Initialize();


/* -------------------  NOT INCLUDING THIS DISTRIBUTED STUFF FOR NOW! -----------------------  

  // Server Routines
  long _compute_smallslices(const unsigned long SLICE_MIN);
  //  bool _read_from_distributed_host(const char* host);
  bool _SliceFileExists(const long slice, const char* host) const;
  bool _incorporate_temp_slice_results(const long slice, const char* host);
  bool _delegate_slice(distributed_info & D, long & slice, long i) const;
  bool _save_current_distributed_state(const distributed_info & D) const;

  bool _SliceFileExists_new(const char* host) const;
  bool _incorporate_temp_slice_results_new(const char* host);
  bool _delegate_slice_new(distributed_info & D, const vector<long> & slice, long i) const;



  // Client Routines
  void _compute_client_slice(const long slice_num); 
  void _WriteTernaryThetaBinary_client(char* hostname, long slice);
  vector<long> _read_slicevector(const char* host);
 
  void _WriteTernaryThetaBinary_client_new(char* hostname);


  // Both
  long _compute_with_slices(const double Cholesky[], const unsigned long SLICE_MIN, const vector<long> slice_vec);

*/


};






// Testing Routines
// ----------------

// Make a set of ternary and quaternary forms to use to test things.



// Testing the Read/WriteTernaryThetaBinary Routines
void Test__BooleanTernaryTheta_read_write();

// Test the Magma read against the compute to check the computations are correct
void Test__MagmaRead_vs_compute();


/* -------------------  NOT INCLUDING THIS DISTRIBUTED STUFF FOR NOW! -----------------------

// Test the compute routine against the compute_distributed routine also
void Test__compute_vs_compute_distributed(); 

// Test the compute routine against the compute_distributed routine also
void Test__compute_vs_multislice_compute_distributed();

*/


// ---------------------------------------

// This has the routine FileExists()
#include "../Utilities/file_utilities.h" 

// These are the boolean_theta_class files
#include "boolean_theta_class__defn.cc"  // Mostly testing routines
#include "boolean_theta_class__methods.cc"

/* This one is not included for now since we don't use it!
#include "boolean_theta_class__compute.cc"
*/

/*
#include "boolean_theta_class__distributed_methods2.cc"
#include "../distributed_theta__client_more.cc"
#include "boolean_theta_class__distributed_methods_new.cc"
#include "../multi-slice_client__more.cc"
*/

#endif
