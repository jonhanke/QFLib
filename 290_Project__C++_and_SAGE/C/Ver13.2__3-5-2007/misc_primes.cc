



// =============================================================================================================
//                                              Use_Next_Primefile(...)
// =============================================================================================================


void Use_Next_Primefile() {
  /*
    Increment the global index primefile_index that says which primefile 
    from primefile_filelist we are currently using, and read in the 
    associated primefile to the vector<long> Big_Prime_List.
  */

  // Declare the external global variables relating to the primes
  extern const char PRIME_DIR[14];
  extern const string PRIMEFILE_LIST[2];
  extern const long NUMBER_OF_PRIMEFILES;
  extern long PRIMEFILE_INDEX;
  extern vector<long> Big_Prime_List;

  
  // SANITY CHECK: Check that we have not exceeded the maximum allowed primefile index
  if (PRIMEFILE_INDEX >= (NUMBER_OF_PRIMEFILES - 1)) {
    cout << "ERROR: We have have already used all primefiles, so we cannot add more!" << endl;
    exit(1);
  }

  // Increment the primefile index
  PRIMEFILE_INDEX += 1;

  
  // Read the shell variable QFLIB_PRIMES_BASE_DIR used to locate the current project output files.
  // (This used to be hardcoded as the const string PRIME_DIR above.)
  const char * PRIMES_BASE_DIR_cstr;
  PRIMES_BASE_DIR_cstr = getenv("QFLIB_PRIMES_BASE_DIR");
  if (PRIMES_BASE_DIR_cstr != NULL) {
    printf ("\n\nThe primes directory shell variable QFLIB_PRIMES_BASE_DIR is: %s\n\n", PRIMES_BASE_DIR_cstr);
    string PRIMES_BASE_DIR(PRIMES_BASE_DIR_cstr);
  } 
  else {
    PRIMES_BASE_DIR_cstr = PRIME_DIR;
  }

  // DIAGNOSTIC
  string PRIMES_BASE_DIR(PRIMES_BASE_DIR_cstr);
  cout << endl << endl;
  cout << "The current primefile directory PRIMES_BASE_DIR is: " << PRIMES_BASE_DIR;
  cout << endl << endl;
  
    
  // Read in a list of primes:                     // Note: These are needed for representability.h
  // -------------------------
  cout << " Starting to read the primes: " << endl;
  cout << " ---------------------------- " << endl;
  PrintTime();
  
  char prime_filename[200]; 
  sprintf(prime_filename, "%s%s", PRIMES_BASE_DIR_cstr, PRIMEFILE_LIST[PRIMEFILE_INDEX].c_str());   
  cout << "  Using prime_filename = " << prime_filename << endl;   
  
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
    
}

