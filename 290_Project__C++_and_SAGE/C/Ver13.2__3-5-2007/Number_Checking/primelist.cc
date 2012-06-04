

 
// Define the Kronecker symbol (a/p)
//  (assuming p is a prime >= 2)
int KroneckerSymbol(long a, long p) {

  // Warning: The modulus symbol retains the sign of the starting number... =|

  if (p==2)
    if (a % 2 == 0)
      return 0;
    else {
      if (abs(a) % 8 == 1)
	return 1;
      else
	return -1;
    }

  else {
    // Ensure a1 is the positive residue of a (mod p)
    long a1 = a % p;
    if (a1 < 0)
      a1 = a1 + p;

    long b = 1;
    for(long i=1;  i<=((p-1)/2); i++)
      b = (b * a1) % p;

    if (b == (p-1))         // Check if b == p-1
      return -1;
    else 
      return int(b);

  }
}



// =======================================================================================================



// Define the function F4primeiso which computes F4(p)
//   (assuming p is prime and isotropic)
double F4PrimeIsoExact(long p, long N, int charvalue) {

  double f4;
  double pp = double(p);  // Warning: This conversion is necessary to prevent a subtle error...(see below)
  f4 = sqrt(pp) * 0.5;
  if ((charvalue == -1) && ((N % p) != 0))
    f4 = f4 * (pp*pp - pp + 1) / (pp*pp + 1);

  return f4;
}




// Define the function F4primeiso which computes F4(p)
//   (assuming p is prime and isotropic)
double F4PrimeIsoApprox(long p, long N, int charvalue) {

  double f4;
  double pp = double(p);  // Warning: This conversion is necessary to prevent a subtle error...(see below)
  f4 = sqrt(pp) * 0.5;
  if ((charvalue == -1) && ((N % p) != 0))
    f4 = f4 - f4/pp;

  return f4;
}



// Define the function F4primeiso which computes F4(p)
//   (assuming p is prime and isotropic)
double F4PrimeIsoUniversal(long p, long N, int charvalue) {

  double f4;
  double pp = double(p);  // Warning: This conversion is necessary to prevent a subtle error...(see below)
  f4 = sqrt(pp) * 0.5;
  f4 = f4 - f4/pp;

  return f4;
}




// =======================================================================================================



// Gives a list of possible primeswith F4(p) <= Bound
//   (Note that we pass the Plist and Flist as arguments here...) 
void PrimeListFromBound(vector<long> & prime_list, vector<double> & F4_value_list, double bound, long N, long chi_top, vector<long> Big_Prime_List) {
  /*
  vector<long> prime_list;
  vector<double> F4_value_list;
  */

  cout << " Entering PrimeListFromBound with: " << endl;
  cout << "       bound = " << bound << endl;
  cout << "       level = " << N << endl;
  cout << "     chi_top = " << chi_top << endl;


  // Create the output filenames: 
  //    primes:  "eligible_primes__bound_N_chitop.data"
  //     F4(p):  "eligible_F4(p)__bound_N_chitop.data"
  char primefilename[200], F4filename[200];    
  sprintf(primefilename, "eligible_primes__b=%f_N=%li_chi=%li__.data", bound, N, chi_top); 
  sprintf(F4filename, "eligible_F4(p)__b=%f_N=%li_chi=%li__.data", bound, N, chi_top); 
  /*
  cout << "\n Filename: " << primefilename << "X" << endl << endl;  
  cout << "\n Filename: " << F4filename << "X" << endl << endl;  
  */


  // Check if these files already exist...
  bool Files_exist = FileExists(primefilename) && FileExists(F4filename);

  // If they exist, then read them in and exit. =)
  if (Files_exist == true) {
    prime_list = ReadVector_long(primefilename);
    F4_value_list = ReadVector_double(F4filename);

    cout << " Read the eligible primes and the F4 list from the datafiles." << endl;

    // If they have the same number of elements then exit, otherwise proceed below...
    if (prime_list.size() == F4_value_list.size())
      return;
    else 
      cout << " PrimeListFromBound Error:  The # of primes and the # of F4(p) don't match! =(" << endl;
  }
  
  // otherwise, continue and regenerate the files...
  else 
    cout << " PrimeListFromBound Error:  Can't open eligible prime and F4(p) files, so we'll regenerate them. =) " << endl;



  // Open the output files and check they're opened correctly
  ofstream primefile, F4file; 
  primefile.open(primefilename);
  F4file.open(F4filename);

  if (! primefile.is_open())
    { cout << "Error opening output file " << primefilename; exit (1); }
  if (! F4file.is_open())
    { cout << "Error opening output file " << F4filename; exit (1); }
  

  // Look for eligible primes
  unsigned long prime_ptr = 0;
  unsigned long p = Big_Prime_List[prime_ptr];      // This sets p=2

  int char_val = KroneckerSymbol(chi_top, p);
  double F4_val = F4PrimeIsoExact(p, N, char_val);

  while ((F4_val <= bound) && (prime_ptr < Big_Prime_List.size() - 1)) {
    prime_list.push_back(p);
    F4_value_list.push_back(F4_val);
    primefile << p << ", ";
    F4file << F4_val << ", ";

    prime_ptr = prime_ptr + 1;
    p = Big_Prime_List[prime_ptr];
                                                                                                                    
    char_val = KroneckerSymbol(chi_top, p);
    F4_val = F4PrimeIsoExact(p, N, char_val);
  }

  cout << " Exiting PrimeListFromBound with: " << endl;
  cout << "           p = " << p << endl;
  cout << "       F4(p) = " << F4_val << endl;
  cout << "   prime_ptr = " << prime_ptr << endl;


  if (prime_ptr >= Big_Prime_List.size() - 1) {
    cout << "ERROR in PrimeListFromBound: We ran out of primes... =(";
    exit(0);
  }                                                                                                        


  // Remove the trailing ", " and insert an EOF
  long file_ptr;
  file_ptr = primefile.tellp();
  primefile.seekp(file_ptr - 2);
  primefile.put(' ');
  primefile.put(EOF);
  primefile.close();

  file_ptr = F4file.tellp();
  F4file.seekp(file_ptr - 2);
  F4file.put(' ');
  F4file.put(EOF);
  F4file.close();

}




// Gives a list of possible primeswith F4(p) <= Bound
//   (Note that we pass the Plist and Flist as arguments here...) 
void PrimeListFromBoundUniversal(vector<long> & prime_list, vector<double> & F4_value_list, double bound, vector<long> Big_Prime_List) {
  /*
  vector<long> prime_list;
  vector<double> F4_value_list;
  */

  cout << " Entering PrimeListFromBound with: " << endl;
  cout << "       bound = " << bound << endl;
  //  cout << "       level = " << N << endl;
  //  cout << "     chi_top = " << chi_top << endl;


  // Create the output filenames: 
  //    primes:  "eligible_primes__bound.data"
  //     F4(p):  "eligible_F4(p)__bound.data"
  char primefilename[60], F4filename[60];  ;  
  sprintf(primefilename, "eligible_primes__b=%f__.data", bound); 
  sprintf(F4filename, "eligible_F4(p)__b=%f__.data", bound); 
  /*
  cout << "\n Filename: " << primefilename << "X" << endl << endl;  
  cout << "\n Filename: " << F4filename << "X" << endl << endl;  
  */


  // Check if these files already exist...
  bool Files_exist = FileExists(primefilename) && FileExists(F4filename);

  // If they exist, then read them in and exit. =)
  if (Files_exist == true) {
    prime_list = ReadVector_long(primefilename);
    F4_value_list = ReadVector_double(F4filename);

    cout << " Read the eligible primes and the F4 list from the datafiles." << endl;

    // If they have the same number of elements then exit, otherwise proceed below...
    if (prime_list.size() == F4_value_list.size())
      return;
    else 
      cout << " PrimeListFromBound Error:  The # of primes and the # of F4(p) don't match! =(" << endl;
  }
  
  // otherwise, continue and regenerate the files...
  else 
    cout << " PrimeListFromBound Error:  Can't open eligible prime and F4(p) files, so we'll regenerate them. =) " << endl;



  // Open the output files and check they're opened correctly
  ofstream primefile, F4file; 
  primefile.open(primefilename);
  F4file.open(F4filename);

  if (! primefile.is_open())
    { cout << "Error opening output file " << primefilename; exit (1); }
  if (! F4file.is_open())
    { cout << "Error opening output file " << F4filename; exit (1); }
  

  // Look for eligible primes
  unsigned long prime_ptr = 0;
  unsigned long p = Big_Prime_List[prime_ptr];      // This sets p=2

  //  int char_val = KroneckerSymbol(chi_top, p);
  double pp = double(p);  // Warning: This conversion is necessary to prevent a subtle overflow error...
  double F4_val = sqrt(pp) * 0.5;
  F4_val = F4_val - F4_val/pp;  // This is the approximate F4 value which is below F4(p). =)

  while ((F4_val <= bound) && (prime_ptr < Big_Prime_List.size() - 1)) {
    prime_list.push_back(p);
    F4_value_list.push_back(F4_val);
    primefile << p << ", ";
    F4file << F4_val << ", ";

    prime_ptr = prime_ptr + 1;
    p = Big_Prime_List[prime_ptr];
                                                                                                                    
    //char_val = KroneckerSymbol(chi_top, p);
    pp = double(p);
    F4_val = sqrt(pp) * 0.5;
    F4_val = F4_val - F4_val/pp;
  }

  cout << " Exiting PrimeListFromBound with: " << endl;
  cout << "           p = " << p << endl;
  cout << "       F4(p) = " << F4_val << endl;
  cout << "   prime_ptr = " << prime_ptr << endl;


  if (prime_ptr >= Big_Prime_List.size() - 1) {
    cout << "ERROR in PrimeListFromBound: We ran out of primes... =(";
    exit(0);
  }                                                                                                        


  // Remove the trailing ", " and insert an EOF
  long file_ptr;
  file_ptr = primefile.tellp();
  primefile.seekp(file_ptr - 2);
  primefile.put(' ');
  primefile.put(EOF);
  primefile.close();

  file_ptr = F4file.tellp();
  F4file.seekp(file_ptr - 2);
  F4file.put(' ');
  F4file.put(EOF);
  F4file.close();

}






// Gives a list of possible primeswith F4(p) <= Bound
//   (Note that we pass the Plist and Flist as arguments here...) 
void PrimeListFromBoundOptimized(vector<long> & prime_list, vector<double> & F4_value_list, double bound, long N, long chi_top, vector<long> Big_Prime_List) {
  /*
  vector<long> prime_list;
  vector<double> F4_value_list;
  */

  unsigned long Fast_F4_transition = 10000;    // Compute the F_4(p) bound exactly for all primes p < 10,000


  cout << " Entering PrimeListFromBound with: " << endl;
  cout << "       bound = " << bound << endl;
  cout << "       level = " << N << endl;
  cout << "     chi_top = " << chi_top << endl;


  // Create the output filenames: 
  //    primes:  "eligible_primes__bound_N_chitop.data2"
  //     F4(p):  "eligible_F4(p)__bound_N_chitop.data2"
  char primefilename[60], F4filename[60];  ;  
  sprintf(primefilename, "eligible_primes__b=%f_N=%ld_chi=%ld__.data2", bound, N, chi_top); 
  sprintf(F4filename, "eligible_F4(p)__b=%f_N=%ld_chi=%ld__.data2", bound, N, chi_top); 
  /*
  cout << "\n Filename: " << primefilename << "X" << endl << endl;  
  cout << "\n Filename: " << F4filename << "X" << endl << endl;  
  */


  // Check if these files already exist...
  bool Files_exist = FileExists(primefilename) && FileExists(F4filename);

  // If they exist, then read them in and exit. =)
  if (Files_exist == true) {
    prime_list = ReadVector_long(primefilename);
    F4_value_list = ReadVector_double(F4filename);

    cout << " Read the eligible primes and the F4 list from the datafiles." << endl;

    // If they have the same number of elements then exit, otherwise proceed below...
    if (prime_list.size() == F4_value_list.size())
      return;
    else 
      cout << " PrimeListFromBound Error:  The # of primes and the # of F4(p) don't match! =(" << endl;
  }
  
  // otherwise, continue and regenerate the files...
  else 
    cout << " PrimeListFromBound Error:  Can't open eligible prime and F4(p) files, so we'll regenerate them. =) " << endl;



  // Open the output files and check they're opened correctly
  ofstream primefile, F4file; 
  primefile.open(primefilename);
  F4file.open(F4filename);

  if (! primefile.is_open())
    { cout << "Error opening output file " << primefilename; exit (1); }
  if (! F4file.is_open())
    { cout << "Error opening output file " << F4filename; exit (1); }
  

  // Look for eligible primes
  unsigned long prime_ptr = 0;
  unsigned long p = Big_Prime_List[prime_ptr];      // This sets p=2

  int char_val = KroneckerSymbol(chi_top, p);
  double F4_val = F4PrimeIsoExact(p, N, char_val);

  // Compute them the slow way for primes < Fast_f4_transitons (= 10,000)
  while ((p <= Fast_F4_transition) && (F4_val <= bound) && (prime_ptr < Big_Prime_List.size() - 1)) {
    prime_list.push_back(p);
    F4_value_list.push_back(F4_val);
    primefile << p << ", ";
    F4file << F4_val << ", ";

    prime_ptr = prime_ptr + 1;
    p = Big_Prime_List[prime_ptr];
                                                                                                                    
    char_val = KroneckerSymbol(chi_top, p);
    F4_val = F4PrimeIsoExact(p, N, char_val);
  }


  // Then compute the fast approximation of F4(p) for the larger primes! =) 
  double pp = double(p);  // Warning: This conversion is necessary to prevent a subtle overflow error...
  F4_val = sqrt(pp) * 0.5;
  F4_val = F4_val - F4_val/pp;  // This is the approximate F4 value which is below F4(p). =)

  while ((F4_val <= bound) && (prime_ptr < Big_Prime_List.size() - 1)) {
    prime_list.push_back(p);
    F4_value_list.push_back(F4_val);
    primefile << p << ", ";
    F4file << F4_val << ", ";

    prime_ptr = prime_ptr + 1;
    p = Big_Prime_List[prime_ptr];
                                                                                                                    
    pp = double(p);
    F4_val = sqrt(pp) * 0.5;
    F4_val = F4_val - F4_val/pp;
  }


  cout << " Exiting PrimeListFromBound with: " << endl;
  cout << "           p = " << p << endl;
  cout << "       F4(p) = " << F4_val << endl;
  cout << "   prime_ptr = " << prime_ptr << endl;


  if (prime_ptr >= Big_Prime_List.size() - 1) {
    cout << "ERROR in PrimeListFromBound: We ran out of primes... =(";
    exit(0);
  }                                                                                                        


  // Remove the trailing ", " and insert an EOF
  long file_ptr;
  file_ptr = primefile.tellp();
  primefile.seekp(file_ptr - 2);
  primefile.put(' ');
  primefile.put(EOF);
  primefile.close();

  file_ptr = F4file.tellp();
  F4file.seekp(file_ptr - 2);
  F4file.put(' ');
  F4file.put(EOF);
  F4file.close();

}


