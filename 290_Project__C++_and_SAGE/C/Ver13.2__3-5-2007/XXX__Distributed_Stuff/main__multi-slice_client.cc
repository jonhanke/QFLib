

#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <gmp.h>
#include <fstream>
#include <math.h>

#include <valarray>


#include <time.h>
#include <string.h>
//#include <stdio.h>




// Additional Libraries from the Older Project...
#include <list>
#include <valarray>
#include <string>

#include <sstream>
#include <iomanip>



using namespace std;


// Note: Since a long has 4 bytes, it can store numbers up to 2^24 = 4,294,967,296
// Our largest prime is less than 15 million, so we only need the mpz_class precision 
// for the square-free numbers. =)




// New files to include

#include <unistd.h>  // Uses gethostname();
#include <stdlib.h>  // Uses atol();






//////////////////////////////////////////
// Global constants used to store files //
//////////////////////////////////////////

// Says where the list of primes are stored
const char PRIME_DIR[] = "/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Primes/";

/*
// Says where the theta functions are
const char THETA_DIR[] = "/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/Theta_Data/";

// Says where the intermediate computations are stored
const char DIST_DIR[] = "/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/TEMP_Distributed_Dir/";

// Says where the ternary_theta_client routine is
//const char EXEC_DIR[] = "/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/Ver7__8-19-04/";
const char EXEC_DIR[] = "/home/postdoc/jonhanke/";
*/



/*
#include "max-min.cc"

#include "LCM.cc"   // This is needed for misc_utilities.cc, but it's repeated in the local densities code!

#include "misc_utilities.h"
#include "misc_utilities.cc"
*/


#include "vec_utils.cc"

/*
#include "cholesky_distributed.cc"

#include "distributed_class.cc"
#include "boolean_theta_class.cc"
#include "boolean_theta_class__methods.cc"
#include "boolean_theta_class__distributed_methods2.cc"
#include "boolean_theta_class__distributed_methods_new.cc"

#include "distributed_theta__client_more.cc"
#include "multi-slice_client__more.cc"
*/

#include "boolean_theta_class.h"







int 
main(int argc, char* argv[])
{
  
  ///////////////////////////////////////////////////////////////////
  // The arguments should be Q[0], ..., Q[5], Precision, slice_num //
  ///////////////////////////////////////////////////////////////////

  // Parse the command line arguments
  cout << " argc = " << argc << endl;
  for(int index=0; index < argc; index++){
    cout << " argv[" << index << "] = " << argv[index] << endl;
  }
  
  // Check the number of arguments
  if (argc != 8) 
    cout << "Error in the number of arguments passed!" << endl;

  // Read in the arguments
  long QQ[6];
  unsigned long Ternary_Precision;
  long slice_number;
  for (long i=0; i < 6; i++)
    QQ[i] = atol(argv[i+1]);
  Ternary_Precision = atol(argv[7]);

  
  // Print out the arguments
  cout << endl << endl;
  for (long i=0; i < 6; i++)
    cout << " QQ[" << i << "] = " << QQ[i] << endl;
  cout << " Ternary_Precision = " << Ternary_Precision << endl;
  cout << endl << endl;
  

  // Find the size of a long
  cout << " sizeof(long) = " << sizeof(long) << endl;



  // Initialize the boolean_ternary_theta class
  boolean_ternary_theta TEMP(QQ, Ternary_Precision);
  
  // Compute the theta function of that slice
  TEMP.compute_client_new();
  


  cout << endl << endl;
  PrintTime();

}
