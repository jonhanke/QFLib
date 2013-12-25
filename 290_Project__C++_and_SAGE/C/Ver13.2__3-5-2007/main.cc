



/////////////////////////
// Basic C/C++ Headers //
/////////////////////////

#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <gmp.h>
#include <fstream>
#include <math.h>
#include <stdlib.h>   // For ultoa() in Matrix_mpz::_GetEisData()
#include <cstdlib>

#include <valarray>


#include <time.h>
#include <string.h>


// Additional Libraries from the Older Project...
#include <list>
#include <valarray>
#include <string>

#include <sstream>
#include <iomanip>

#include <climits>



////////////////
// Namespaces //
////////////////

using namespace std;





//////////////////////////////////////////////////////////////////
// Set Global variables related to maintaining a list of primes //
//////////////////////////////////////////////////////////////////

//const char PRIME_DIR[14] = "../../Primes/";   // This gives the default directory for where the primefiles are stored,

const string PRIME_DIR_STRING = "../../Primes/";   // This gives the default directory for where the primefiles are stored,
                                                  // which is overridden by the shell variable QFLIB_PRIME_DIR.

// Define a list of space-delimited text files containing the primes up to some bound 
// in increasing order, with the list of files in order of increasing size.
const string PRIMEFILE_LIST[3] =
  {
    string("primes_lt_1000000.txt"),
    string("primelist_upto_15millionth_prime.txt"),
    string("primelist_upto_40millionth_prime.txt")
  };
const long NUMBER_OF_PRIMEFILES = 3;

// Define a global index pointing to the current primefile in use from primefile_filelist
long PRIMEFILE_INDEX = -1;

  // Note: Since a long has 4 bytes, it can store numbers up to 2^24 = 4,294,967,296
  // Our largest prime is less than 15 million, so we only need the mpz_class precision 
  // for the square-free numbers. =)


// The global list of primes used throughout.
vector<long> Big_Prime_List;






///////////////////
// My QF headers //
///////////////////

// Put the headers here!
#include "GMP_class_extras/mpz_class_extras.h"
#include "Matrix_mpz/Matrix_mpz.h"

#include "New_Boolean_Theta_Class/boolean_theta_class.h"  



#include "Matrix_mpz/Theta.h"   // WARNING: The "Theta.cc" file is *not* included in this!


/*
// Local Diagnostic Routines -- to be moved into the Matrix_mpz class...     
#include "Siegel_Diagnostic/siegel_diagnostic.h"
*/


/*
#include "Siegel_Diagnostic/readseries.cc"        // TO DO: This should be replaced by a theta series class...
*/




// #include "misc_utilities.h"   // <<--- Do we need this?

// #include "distributed_class.cc"  // <<--- What's this?



#include "Utilities/vec_utils.h"
#include "Utilities/file_utilities.h"

#include "Number_Checking/primelist.h"
#include "Number_Checking/checknums.h"
#include "Number_Checking/squarefreelist.h"

#include "local_condition_class.h"




#include "scraps5.h"  // Tests the modified densities on imprimitive forms!



#include "qf_project_class.h"
#include "qf_datafiles_class.h"
#include "representability.h"  



// To Do: Look at this next...
// #include "new_scraps.h"  


#include "power_series.h"
#include "local_checking_stats.h"
#include "misc_primes.h"



// PARI Stuff:
// -----------
#include <pari.h>    // BEWARE: Namespace conflicts!!!
#include "Matrix_mpz/PARI.cc"   // This has the conversion routine! =)
#include "Matrix_mpz/Theta.cc"   // WARNING: The "Theta.h" file is included above, and does *not* include this!





 
// =============================================================================================================
//                                              main(...)
// =============================================================================================================


int main(int argc, char *argv[])
{

  // /*
  // Testing the size of certain data-types
  long test_long;
  unsigned long test_unsigned_long;
  double test_double;
  //  long double test_long_double;              // WARNING: SOMETIMES THE 2nd LINE BELOW MAKES A SEGFAULT!!!

  cout << " The size of LONG is: " << sizeof(test_long) << endl;
  cout << " The size of UNSIGNED LONG is: " << sizeof(test_unsigned_long) << endl;
  cout << " The size of DOUBLE is: " << sizeof(test_double) << endl;
  //  cout << " The size of LONG DOUBLE is: " << sizeof(test_long_double) << endl;
  cout << endl;

  cout << " The default LONG is: " << test_long << endl;
  cout << " The default UNSIGNED LONG is: " << test_unsigned_long << endl;
  cout << " The default DOUBLE is: " << test_double << endl;
  // cout << " The default LONG DOUBLE is: " << test_long_double << endl;
  cout << endl;

  cout << " The smallest LONG is: " << LONG_MIN << endl;
  cout << " The largest LONG is: " << LONG_MAX << endl;
  cout << " The largest UNSIGNED LONG is: " << ULONG_MAX << endl;
  //  cout << " The largest DOUBLE is: " << DOUBLE_MAX << endl;
  //  cout << " The largest LONG DOUBLE is: " << test_long_double << endl;
  cout << endl;

  vector<unsigned long> test_vec_UL;
  cout << " The largest a vector<unsigned long> can be is: " << test_vec_UL.max_size() << endl;
  cout << endl;
  // */



  // Initialize the Pari Stack with 4MB and all primes up to 2
  pari_init(4000000,2);

  
  // Declare some global filenames
  //  char* prime_file; 
  char* theta_file;
  
  
  
  // =================================================================================================================
  // =================================================================================================================


  
  
  // Timestamp for program start:
  // ----------------------------
  cout << " ================================================" << endl;
  cout << " Program started at: ";
  PrintTime();
  cout << " ================================================" << endl << endl;


  // Read the command line arguments:  
  //      0 = "./main"??
  //      1 = project_name 
  //      2 = prime_list     // Depricated -- no longer used!
  //      3 = form_file
  //      4 = Cusp_Const_Dir
  //      5 = cusp_const_prefix
  //      6 = form_number
  //      7 = form_number (optional -- last one in a range if it's included)
  cout << " We read " << argc << " command line arguments! ";
  assert((argc == 7) || (argc == 8));
    
  char * pEnd;
  string CL_project_name(argv[1]);
  string CL_prime_filename(argv[2]);   // Depricated -- no longer used!
  string CL_form_filename(argv[3]);
  string CL_cusp_const_dir(argv[4]);
  string CL_cusp_const_prefix(argv[5]);
  long form_number = strtol(argv[6], &pEnd, 10);
  bool form_range_flag = false;
  long last_form_number;
  if (argc == 8) {
    form_range_flag = true;
    last_form_number = strtol(argv[7], &pEnd, 10);
  }

  assert(form_number >= 1);
  /*
  assert(form_number <= 6560);
  */

 

    
  // DIAGNOSTIC  
  cout << "   Read in CL_project_name as: " << CL_project_name << endl;
  cout << " Read in (but not using) CL_prime_filename as: " << CL_prime_filename << endl;
  cout << "           Read the form # as: " << form_number << endl;
  cout << endl;



  // ================================== Reading in the primes =======================================


  // Load in the initial primefile
  //cout << "Before Use_Next_Primefile()" << endl;
  Use_Next_Primefile();
  //cout << "After Use_Next_Primefile()" << endl;
  
  // Note: Since a long has 4 bytes, it can store numbers up to 2^24 = 4,294,967,296
  // Our largest prime is less than 15 million, so we only need the mpz_class precision 
  // for the square-free numbers. =)



  // =============================== Find/create the project data directory ===============================

  // Set the external data directory variable to use throughout the program
  string ABSOLUTE_PROJECT_PATH_STRING = GetAbsoluteProjectDirPath();



  // ======================================= Exception Finding ==============================================


  //  string form_file("~/290_Project/290_Cusp_info/data2.m");
  //  string form_file = CL_form_filename;



  // Find exceptions using the QF_Project class
  string command_line_project_name(CL_project_name);
  /*
  string project_name_8_15_2005("290_Project__8-15-2005__with_Approx_thetas");
  string project_name_8_11_2005("290_Project__8-11-2005__with_Approx_thetas");
  string project_name_8_8_2005("290_Project__8-8-2005__with_Approx_thetas");
  string project_name_7_14_2005("290_Project__7-14-2005__with_Approx_thetas");
  string project_name_5_16_2005("290_Project__5-16-2005__with_Approx_thetas");
  string project_name_Final("290_Project__Final");
  string project_name2("290_Project__2");
  string project_name("290_Project");
  */



  // Define some excluded forms for this project
  set<long> Excluded_Forms;

  /*
  Excluded_Forms.insert(5875);                // Big split local cover (D=697) 
  Excluded_Forms.insert(6414);                // This is the big one! =)
  */

  /*
  Excluded_Forms.insert(258);                 // Not Locally universal! =(    
  Excluded_Forms.insert(757);                 // Not Locally universal! =(    
  Excluded_Forms.insert(1825);                // Not Locally universal! =(    
  Excluded_Forms.insert(3536);                // Not Locally universal! =(     (See above for details...)
  Excluded_Forms.insert(3819);                // Not Locally universal! =(    
  */


  // Define the project
  QF_Project My_Project(command_line_project_name, CL_form_filename, CL_cusp_const_dir, CL_cusp_const_prefix, Excluded_Forms);


  /*
  // Set the bounds to check 
  long first = 1;
  long last = 6560;
  
  // Find all exceptions (for non-excluded forms) within those bounds
  cout << endl << endl;
  cout << " About to check Form #" << first << " to Form #" << last << endl;
  cout << endl << endl; 
  */



  // Compute the exceptions and print the results
  if (form_range_flag == true) {
    My_Project.FindExceptions("minimal", form_number, last_form_number);
    My_Project.PrintExceptions(form_number, last_form_number);
  }
  else {
    My_Project.FindExceptions("minimal", form_number, form_number);
    My_Project.PrintExceptions(form_number, form_number);
  }




  /* 
  // Testing the uniform checking method! =) 
  string project_name__uniform1("290_Project__uniform1");
  QF_Project My_Project(project_name__uniform1, form_file);
  My_Project.FindExceptions("uniform");
  */




  // Timestamp for program end:
  // --------------------------
  cout << " ================================================" << endl;
  cout << " Program finished at: ";
  PrintTime();
  cout << " ================================================" << endl << endl;


}
