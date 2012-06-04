



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

#include <valarray>


#include <time.h>
#include <string.h>
//#include <stdio.h>
//#include <stdlib.h>


// Additional Libraries from the Older Project...
#include <list>
#include <valarray>
#include <string>

#include <sstream>
#include <iomanip>





////////////////
// Namespaces //
////////////////

using namespace std;

namespace JON { }
using namespace JON;




//////////////////////////////////////////
// Global constants used to store files //
//////////////////////////////////////////

// Says where the list of primes are stored
const char PRIME_DIR[] = "/home/postdoc/jonhanke/290_Project/Primes/";

// Says where the project is stored
//const char ABSOLUTE_PROJECT_PATH[] = "/ytmp/QF_Project_Data/";   // This is for Austin and grid
extern const char ABSOLUTE_PROJECT_PATH[23] = "/ztmp/QF_Project_Data/";   // This is for the gridX machines







///////////////////
// My QF headers //
///////////////////

// Put the headers here!
#include "GMP_class_extras/mpz_class_extras.h"
#include "Matrix_mpz/Matrix_mpz.h"

//#include "Boolean_Ternary_Theta_Class/boolean_ternary_theta_class.h"    // <------- This is on it's way out!

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


/*
// These are "ad hoc" routines which check representability of ternary forms.
// (They need to be updated and consolidated!)
#include "Ternary_Stuff/local_representability.h"
#include "Ternary_Stuff/new_local_representability.h"
#include "Ternary_Stuff/scraps3.h"
*/


#include "local_condition_class.h"





/*
#include "Number_Checking/checking.h"
#include "scraps4.h"
*/

#include "scraps5.h"  // Tests the modified densities on imprimitive forms!



#include "qf_project_class.h"
#include "qf_datafiles_class.h" // I think this last one is unnecessary now, since it's included (and needs) QF_Project
#include "representability.h"   // I think this last one is unnecessary now, since it's included (and needs) QF_Project



// To Do: Look at this next...
// #include "new_scraps.h"  



#include "power_series.h"

#include "local_checking_stats.h"







/*
// Put the real routines here
//#include "misc_utilities.cc"
#include "vec_utils.cc"
#include "file_stuff.cc"
#include "primelist.cc"
#include "checknums.cc"
#include "squarefreelist.cc"
#include "readmagmaseries.cc"
//#include "maketheta.cc"
#include "maketheta_small.cc"

//#include "more_scraps.cc"  --  What happened to this? =|


//#include "anisoprimelist.cc"  // This includes all of the old local density libraries


#include "local_representability.cc"
#include "new_local_representability.cc"
//#include "scraps3.cc"
#include "checking.cc"

#include "local_condition_class.cc"

#include "more_local.cc"
#include "more_local__tests.cc"

#include "scraps4.cc"
*/


/*
#include "boolean_theta_class__methods.cc"
#include "boolean_theta_class__distributed_methods2.cc"
#include "boolean_theta_class__distributed_methods_new.cc"
#include "distributed_theta__client_more.cc"
*/

//#include "boolean_theta_class.h"


//#include "Misc_Tests.cc"



// Big list of primes
vector<long> Big_Prime_List;




/*
// These are still in development
#include "ternary_qf_class.cc"          // <------ Are we even using this anymore??? =|   (moved to folder XXX_UNUSED)
*/

/*  THESE ARE NOT USED FOR NOW...  
#include "Ternary_Exceptions/Ternary_Exceptions.h"           // (moved to folder XXX_UNUSED)
#include "Quaternary_Exceptions/Quaternary_Exceptions.h"     // This is included in representability.cc!
*/



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
  //      2 = prime_list
  //      3 = form_file
  //      4 = Cusp_Const_Dir
  //      5 = cusp_const_prefix
  //      6 = form_number
  //      7 = form_number (optional -- last one in a range if it's included)
  cout << " We read " << argc << " command line arguments! ";
  assert((argc == 7) || (argc == 8));
    
  char * pEnd;
  string CL_project_name(argv[1]);
  string CL_prime_filename(argv[2]);
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
  cout << " Read in CL_prime_filename as: " << CL_prime_filename << endl;
  cout << "           Read the form # as: " << form_number << endl;
  cout << endl;



  // =================================================================================================================




  /*
  Matrix_mpz M3(3,3);
  M3(1,1) = 1;
  M3(1,2) = 2;
  M3(1,3) = 3;
  M3(2,1) = 4;
  M3(2,2) = 5;
  M3(2,3) = 6;
  M3(3,1) = 7;
  M3(3,2) = 8;
  M3(3,3) = 9;

  Matrix_mpz N3(4,4);
  N3 = PARI__to__Matrix_mpz(Matrix_mpz__to__PARI(M3));
  cout << N3 << endl;

  exit(1);
  */




  // ================================== Reading in the primes =======================================
  

  
  // Set the file with all of the primes:
  // ------------------------------------
  //prime_file = "primes_lt_1000000.txt";
  //prime_file = "primelist_upto_15millionth_prime.txt";
  //prime_file = "primelist_upto_40millionth_prime.txt";
  string prime_file = CL_prime_filename;   // Get the filename from the command line! =)



  // Note: Since a long has 4 bytes, it can store numbers up to 2^24 = 4,294,967,296
  // Our largest prime is less than 15 million, so we only need the mpz_class precision 
  // for the square-free numbers. =)



  
  // Read in a list of primes:                                    // Note: These are needed for representability.h
  // -------------------------
  cout << " Starting to read the primes: " << endl;
  cout << " ---------------------------- " << endl;
  PrintTime();


  char prime_filename[200]; 
  sprintf(prime_filename, "%s%s", PRIME_DIR, prime_file.c_str());   


  /*
  // DIAGNOSTIC
  cout << "  prime_file = " << prime_file << endl; 
  cout << "  prime_filename = " << prime_filename << endl; 
  */


  Big_Prime_List = ReadVector_long(prime_filename);
  cout << " Read in " << Big_Prime_List.size() << " prime numbers." << endl;


  // OUTPUT 
  PrintHeadV(Big_Prime_List, 5);
  PrintTailV(Big_Prime_List, 5);
  
  cout << " Finished reading in the primes " << endl << endl;
  PrintTime();



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
