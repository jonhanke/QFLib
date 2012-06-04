

////////////////////////////////////////////////////////////////////
// Compute an upper bound for the precision needed in the ternary //
// theta series to check the eligible square-free numbers.        //
////////////////////////////////////////////////////////////////////

double ThetaPrecision(double B, long N, long chi_top, long diag_coeff) {

  // Entry conditions:
  /*
  cout << " Using B = " << B 
       << "  N = " << N 
       << "  chi_top = " << chi_top 
       << "  diag_coeff = " << diag_coeff << endl;  
  */


  // Deal with the bounds of zero gracefully
  if (B==0)
    return 0;


  //  valarray<long> prime_list;
  const long prime_list[] = 
    {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 
     101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 
     197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 
     311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 
     431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 
     557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 
     661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 
     809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 
     937, 941, 947, 953, 967, 971, 977, 983, 991, 997}; 
  
  const long Prime_List_Length = 168;


  // Find the allowed number of prime factors (up to 100)
  double F4_total = 1;
  long p, i=0;
  while ((F4_total < B) && (i < Prime_List_Length)) {
    p = prime_list[i];    
    F4_total = F4_total * F4PrimeIsoExact(p, N, KroneckerSymbol(chi_top, p)); 
    i++;
  }
  
  // Error Checking:  Check we haven't used all of the primes. =)
  if (i == Prime_List_Length) {
    cout << "ERROR1 in ThetaPrecision:  We've used up all of the primes < 1000. " << endl;
    cout << "  Your bound must have been ridiculously huge! " << endl;
    exit(1);
  }

  long Allowed_primes = i-1;
  

  // Find the product of the first N adjustment factors (N = Max allowed number of primes) 
  long ct = 0;
  double adjustment = 1;
  i = 0; p = prime_list[i];
  while ((ct < Allowed_primes) && (i < Prime_List_Length)) {        // Use < or <= ???
    if (((N % p) != 0) && (KroneckerSymbol(chi_top, p) == -1)) {
      adjustment = adjustment * (1 - 1/double(p+1));
      ct++;
    }
    
    i++; p = prime_list[i];      
  }


  // Error Checking:  Check we haven't used all of the primes or chi_top is 1. =)
  if ((i == Prime_List_Length) && (chi_top != 1)){
    cout << "ERROR2 in ThetaPrecision:  We've used up all of the primes < 1000. " << endl;
    cout << "  Your bound must have been ridiculously huge, or the primes are biased! " << endl;
    exit(1);
  }


  // Compute the needed precision 
  double T_bound;
  T_bound = (pow(2.0, 2.0 * Allowed_primes) * B * B) / (adjustment * adjustment);

  double Diff_bound;
  Diff_bound = 2 * sqrt(T_bound * diag_coeff) + diag_coeff;

  /*
    cout << " 1/adjustment is " << (1/adjustment) << endl;
    cout << " T_bound is " << T_bound << endl;
    cout << " Diff_bound is " << Diff_bound << endl;
  */


  // INSERT A FACTOR OF 5 TO GET THE CORRECT PRECISION FOR A DEPTH 5 HUNT! -- See 8/4/04 notes.
  Diff_bound = 5 * Diff_bound;


  return Diff_bound;

}





/////////////////////////////////////////////////////////////////
// Compute the Choelesky Decomposition of a 3x3 quadratic form //
/////////////////////////////////////////////////////////////////

//ChoeleskyDecomp(mpz_class Q)




////////////////////////////////////
// Find all vectors with Q(x) < C //
////////////////////////////////////

void Short123(vector<mpz_class> & theta, mpz_class C) { 


  cout << endl << "Entering Short123" << endl;
  cout << " Using C = " << C << endl;


  mpz_class CC(C+1);  // This is just C+1, which gives us precision up to C

  // ERROR CHECKING: Check that C+1 fits in an unsigned long.
  if (CC.fits_uint_p() == false) {
    cout << " Error in Short123: The value C+1 = " << CC << " is too big!" << endl;
    exit(1);
  }


  // Resize the theta vector
  theta.resize(0);
  theta.resize(CC.get_ui());


  const mpz_class Q[4][4] = {{0,0,0,0}, {0, 1, 0, 0}, {0, 0, 3, 0}, {0, 0, 0, 5}};  // WARNING: We only use indices 1--3.
  const long n = 3;

  /*
  for (long i=1; i<=n; i++) {
    for (long j=1; j<=n; j++)
      cout << Q[i][j] << " ";
    cout << endl;
  }
  */


  // 1. Initialize
  long i = n;
  vector<mpz_class> T(n+1, 0);  // Note: We use n+1 so we can index the entries as 1 --> n
  vector<mpz_class> U(n+1, 0);
  T[i] = C;
  U[i] = 0;
  double Z;
  vector<mpz_class> L(n+1, 0);
  vector<mpz_class> x(n+1, 0);

  cout << "Here 1" << endl;

  // 2. Compute bounds
  Z = sqrt(T[i].get_d() / Q[i][i].get_d());
  L[i] = mpz_class(floor(Z - U[i].get_d()));  // Check this is ok...
  /*
  cout << " L[i] float gives :     " << floor(Z - U[i].get_d()) << endl;
  cout << " L[i] mpz_class gives : " << L[i] << endl;
  */

  x[i] = mpz_class(ceil(-Z - U[i].get_d()) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i].get_d() -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */
  
  
  bool done_flag = false;
  mpz_class Q_val;

  
  // Big loop which runs through all vectors
  while (done_flag == false) {

    // Loop through until we get to i=1 (so we defined a vector x)
    do {
      
      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }
      
      // 3b. Main loop
      if (i>1) {
	T[i-1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i]);
	i = i - 1;
	for(long j=i+1; j<=n; j++)
	  U[i] = U[i] + Q[i][j] * x[j];
	
	// Now go back and compute the bounds...
	// 2. Compute bounds
	Z = sqrt(T[i].get_d() / Q[i][i].get_d());
	L[i] = mpz_class(floor(Z - U[i].get_d()));  
	x[i] = mpz_class(ceil(-Z - U[i].get_d()) - 1);  
      }
      
    } 
    while (i > 1);
    
    
    // 4. Solution found (This happens when i=1)
    /*
    cout << " x = [ " << x[1] << ", " << x[2] << ", " << x[3] << " ]" << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */
    Q_val = C - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    if (Q_val <= C)
      theta[Q_val.get_ui()] = theta[Q_val.get_ui()] + 2;


   
 
    // 5. Check if x = 0, for exit condition. =)
    long j=1;
    done_flag = true;
    while (j<=n) {
      if (x[j] != 0)
	done_flag = false;
      j++;
    }    
  }

  // Subtract off the doubled count for zero. =)
  theta[0] = theta[0] - 1;
   
  cout << "Leaving Short123" << endl << endl;
  
}







//////////////////////////////////////
// Find all vectors with Q(x) < C   //
// (uses long instead of mpz_class) //
//////////////////////////////////////

void FastShort123(vector<mpz_class> & theta, unsigned long C) { 

  
  cout << endl << "Entering FastShort123" << endl;
  cout << " Using C = " << C << endl;


  // ERROR CHECKING: Check that C+1 fits in an unsigned long.
  if (C+1 < C) {
    cout << " Error in Short123: The value C = " << C << " is too big!" << endl;
    exit(1);
  }


  // Resize the theta vector
  theta.resize(0);
  theta.resize(C+1);


  const long Q[4][4] = {{0,0,0,0}, {0, 1, 0, 0}, {0, 0, 3, 0}, {0, 0, 0, 5}};  // WARNING: We only use indices 1--3.
  const long n = 3;

  /*
  for (long i=1; i<=n; i++) {
    for (long j=1; j<=n; j++)
      cout << Q[i][j] << " ";
    cout << endl;
  }
  */


  // 1. Initialize
  long i = n;
  vector<long> T(n+1, 0);  // Note: We use n+1 so we can index the entries as 1 --> n
  vector<long> U(n+1, 0);
  T[i] = C;
  U[i] = 0;
  double Z;
  vector<long> L(n+1, 0);
  vector<long> x(n+1, 0);

  cout << "Here 1" << endl;

  // 2. Compute bounds
  Z = sqrt(double(T[i]) / double(Q[i][i]));
  L[i] = long(floor(Z - double(U[i])));  // Check this is ok...
  /*
  cout << " L[i] float gives :     " << floor(Z - U[i].get_d()) << endl;
  cout << " L[i] mpz_class gives : " << L[i] << endl;
  */

  x[i] = long(ceil(-Z - double(U[i])) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i].get_d() -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */
  
  
  bool done_flag = false;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...

  
  // Big loop which runs through all vectors
  while (done_flag == false) {

    // Loop through until we get to i=1 (so we defined a vector x)
    do {
      
      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }
      
      // 3b. Main loop
      if (i>1) {
	T[i-1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i]);
	i = i - 1;
	for(long j=i+1; j<=n; j++)
	  U[i] = U[i] + Q[i][j] * x[j];
	
	// Now go back and compute the bounds...
	// 2. Compute bounds
	Z = sqrt(double(T[i]) / double(Q[i][i]));
	L[i] = long(floor(Z - double(U[i])));  
	x[i] = long(ceil(-Z - double(U[i])) - 1);  
      }
      
    } 
    while (i > 1);
    
    
    // 4. Solution found (This happens when i=1)
    /*
    cout << " x = [ " << x[1] << ", " << x[2] << ", " << x[3] << " ]" << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */
    Q_val = C - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    if (Q_val <= C)
      theta[Q_val] = theta[Q_val] + 2;


   
 
    // 5. Check if x = 0, for exit condition. =)
    long j=1;
    done_flag = true;
    while (j<=n) {
      if (x[j] != 0)
	done_flag = false;
      j++;
    }    
  }

  // Subtract off the doubled count for zero. =)
  theta[0] = theta[0] - 1;
   
  cout << "Leaving FastShort123" << endl << endl;
  
}







////////////////////////////////////
// Find all vectors with Q(x) < C //
////////////////////////////////////

void FastShort123Binary(unsigned long theta[], unsigned long C) { 


  // *** WARNING:  We assume that theta has size ceil(C / 32), so it can hold the precision. 

  

  cout << endl << "Entering FastShort123Binary" << endl;
  cout << " Using C = " << C << endl;


  // ERROR CHECKING: Check that C+1 fits in an unsigned long.
  if (C+1 < C) {
    cout << " Error in Short123: The value C = " << C << " is too big!" << endl;
    exit(1);
  }


  // Zero out the theta function
  for (unsigned long i=0; i < (unsigned long) ceil((double) C / (double) 32); i++)
    theta[i] = 0;


  /*
  unsigned long theta[(C >> 5) + 1];   // Allocate at least ceil(C/32) numbers
  */


  const long Q[4][4] = {{0,0,0,0}, {0, 1, 0, 0}, {0, 0, 3, 0}, {0, 0, 0, 5}};  // WARNING: We only use indices 1--3.
  const long n = 3;

  /*
  for (long i=1; i<=n; i++) {
    for (long j=1; j<=n; j++)
      cout << Q[i][j] << " ";
    cout << endl;
  }
  */


  // 1. Initialize
  long i = n;
  vector<long> T(n+1, 0);  // Note: We use n+1 so we can index the entries as 1 --> n
  vector<long> U(n+1, 0);
  T[i] = C;
  U[i] = 0;
  double Z;
  vector<long> L(n+1, 0);
  vector<long> x(n+1, 0);

  cout << "Here 1" << endl;

  // 2. Compute bounds
  Z = sqrt(double(T[i]) / double(Q[i][i]));
  L[i] = long(floor(Z - double(U[i])));  // Check this is ok...
  /*
  cout << " L[i] float gives :     " << floor(Z - U[i].get_d()) << endl;
  cout << " L[i] mpz_class gives : " << L[i] << endl;
  */

  x[i] = long(ceil(-Z - double(U[i])) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i].get_d() -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */
  
  
  bool done_flag = false;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...

  
  // Big loop which runs through all vectors
  while (done_flag == false) {

    // Loop through until we get to i=1 (so we defined a vector x)
    do {
      
      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }
      
      // 3b. Main loop
      if (i>1) {
	T[i-1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i]);
	i = i - 1;
	U[i] = 0;
	for(long j=i+1; j<=n; j++)
	  U[i] = U[i] + Q[i][j] * x[j];
	
	// Now go back and compute the bounds...
	// 2. Compute bounds
	Z = sqrt(double(T[i]) / double(Q[i][i]));
	L[i] = long(floor(Z - double(U[i])));  
	x[i] = long(ceil(-Z - double(U[i])) - 1);  
      }
      
    } 
    while (i > 1);
    
    
    // 4. Solution found (This happens when i=1)
    /*
    cout << " x = [ " << x[1] << ", " << x[2] << ", " << x[3] << " ]" << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */
    Q_val = C - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    if (Q_val <= C) {
      theta[Q_val >> 5] = theta[Q_val >> 5] | (1 << (Q_val % 32));  // Set the appropriate bit if it wasn't set
	//      theta[Q_val] = theta[Q_val] + 2;
    }

   
 
    // 5. Check if x = 0, for exit condition. =)
    long j=1;
    done_flag = true;
    while (j<=n) {
      if (x[j] != 0)
	done_flag = false;
      j++;
    }    
  }

  // Subtract off the doubled count for zero. =)
  //  theta[0] = theta[0] - 1;
   
  cout << "Leaving FastShort123Binary" << endl << endl;


  //  theta_bin = theta;

  //  theta_bin[0] = 12;   // Testing the passing of the array -- it's good! =0
  
}






/////////////////////////////////////////
// Rounds a double to an unsigned long //
// (used below in FastBinaryTheta)     //
/////////////////////////////////////////
#define round(x) ((x)>=0?(unsigned long)((x)+0.5):(unsigned long)((x)-0.5))


////////////////////////////////////
// Find all vectors with Q(x) <= C //
////////////////////////////////////

//! \deprecated  This has been incorporated into boolean_ternary_theta::compute().

void FastBinaryTheta(unsigned long theta[], double QQ[], double C) { 


  // *** WARNING:  We assume that theta has size ceil(C / 32), so it can hold the precision. 

  

  cout << endl << "Entering FastBinaryTheta" << endl;
  cout << " Using C = " << C << endl;


  // ERROR CHECKING: Check that C+1 fits in an unsigned long.


  /*
  unsigned long theta[(C >> 5) + 1];   // Allocate at least ceil(C/32) numbers
  */


  // Zero out the theta function
  for (unsigned long i=0; i < (unsigned long) ceil(C/32); i++)
    theta[i] = 0;


  const long n = 3;

  // Make the (constant) matrix Q from QQ  -- NOTE: We only use indices 1--3.
  double Q[n+1][n+1];
  for (long i=1; i<=n; i++) 
    for (long j=1; j<=n; j++) 
      Q[i][j] = 0;   // Clear the matrix
  long counter = 0;
  for (long i=1; i<=n; i++) 
    for (long j=i; j<=n; j++) {
      Q[i][j] = QQ[counter];  // Put QQ in the upper triangular part
      counter++;
    }



  // Print Q
  //  /*
  cout << "Using Q = " << endl;
  for (long i=1; i<=n; i++) {
    for (long j=1; j<=n; j++)
      cout << Q[i][j] << " ";
    cout << endl;
  }
  //  */


  // 1. Initialize
  long i = n;
  vector<double> T(n+1, 0);  // Note: We use n+1 so we can index the entries as 1 --> n
  vector<double> U(n+1, 0);
  T[i] = C;
  U[i] = 0;
  double Z;
  vector<long> L(n+1, 0);
  vector<long> x(n+1, 0);


  // 2. Compute bounds
  Z = sqrt(T[i] / Q[i][i]);
  L[i] = long(floor(Z - U[i]));  // Check this is ok...
  /*
  cout << " L[i] float gives :     " << floor(Z - U[i].get_d()) << endl;
  cout << " L[i] mpz_class gives : " << L[i] << endl;
  */

  x[i] = long(ceil(-Z - U[i]) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i].get_d() -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */
  
  
  bool done_flag = false;
  double Q_val_double;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...

  
  // Big loop which runs through all vectors
  while (done_flag == false) {

    // Loop through until we get to i=1 (so we defined a vector x)
    do {
      
      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }
      
      // 3b. Main loop
      if (i>1) {
	/*
	cout << " i = " << i << endl;
	cout << " T[i] = " << T[i] << endl;
	cout << " Q[i][i] = " << Q[i][i] << endl;
	cout << " x[i] = " << x[i] << endl;
	cout << " U[i] = " << U[i] << endl;
	cout << " x[i] + U[i] = " << (x[i] + U[i]) << endl;	
	cout << " T[i-1] = " << T[i-1] << endl;		
	*/
	T[i-1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i]);
	/*
	cout << " T[i-1] = " << T[i-1] << endl;		
	cout << endl;
	*/
	i = i - 1;
	U[i] = 0;
	for(long j=i+1; j<=n; j++)
	  U[i] = U[i] + Q[i][j] * x[j];
	
	// Now go back and compute the bounds...
	// 2. Compute bounds
	Z = sqrt(T[i] / Q[i][i]);
	L[i] = long(floor(Z - U[i]));  
	x[i] = long(ceil(-Z - U[i]) - 1);  
      }
      
    } while (i > 1);
    
    
    // 4. Solution found (This happens when i=1)
    /*
    cout << " x = [ " << x[1] << ", " << x[2] << ", " << x[3] << " ]" << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */
    Q_val_double = C - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    Q_val = (unsigned long) round(Q_val_double);

    //    cout << " Float = " << Q_val_double << "   Long = " << Q_val << "  XX " << endl;
    /*
    cout << " The float value is " << Q_val_double << endl;
    cout << " The associated long value is " << Q_val << endl;
    cout << endl;
    */

    if (Q_val <= C) {
      theta[Q_val >> 5] = theta[Q_val >> 5] | (1 << (Q_val % 32));  // Set the appropriate bit if it wasn't set
	//      theta[Q_val] = theta[Q_val] + 2;
    }


 
    // 5. Check if x = 0, for exit condition. =)
    long j=1;
    done_flag = true;
    while (j<=n) {
      if (x[j] != 0)
	done_flag = false;
      j++;
    }    
  }

   
  cout << "Leaving FastBinaryTheta" << endl << endl;

}




////////////////////////////////////////////////////////////////
// Writes the boolean ternary theta function as a binary file //
////////////////////////////////////////////////////////////////

//! \deprecated Moved to the boolean_ternary_theta::write()
void WriteTernaryThetaBinary(unsigned long theta_bin[], long form[6], unsigned long precision) {

  // *** WARNING:  We assume that theta_bin has size precision/32 
  // ***           (which could also be coded as precision >> 5).

  // Notation:
  // ---------
  // precision of the theta function ==> ... + O(x^(precision + 1))
  //
  // Q = [ a, b, c ]
  //     [ 0, d, e ]
  //     [ 0, 0, f ]

  extern char THETA_DIR[];  // This is the global directory for the theta function files
  
  // Make the filename "ternary_bool_theta__a_b_c_d_e_f__precision.bindata"
  char filename[200];  
  sprintf(filename, "%sternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  THETA_DIR, form[0], form[1], form[2], form[3], form[4], form[5], precision); 
  
  // Compute the length of the array of long:  n = ceil(precision / 32)
  unsigned long n;
  n = precision >> 5;
  if ((precision % 32) != 0)
    n = n + 1;
    
  // -------------------------------------------------------------------

  // Open the file for writing and check it opened correctly
  ofstream fileout;
  fileout.open(filename, ios::out | ios::trunc | ios::binary);

  if (! fileout.is_open())
    { cout << "Error opening output file " << filename; exit (1); }

  // Write the array
  fileout.write( (char*) theta_bin, sizeof(long) * n);

  // Close the file
  fileout.close();

}






/////////////////////////////////////////////////////////////////
// Reads the boolean ternary theta function from a binary file //
/////////////////////////////////////////////////////////////////

//! \deprecated Moved to the boolean_ternary_theta::read()
void ReadTernaryThetaBinary(unsigned long theta_read[], long form[6], unsigned long precision) {

  // *** WARNING:  We assume that theta_bin has size precision/32 
  // ***           (which could also be coded as precision >> 5).

  // Notation:
  // ---------
  // precision of the theta function ==> ... + O(x^(precision + 1))
  //
  // Q = [ a, b, c ]
  //     [ 0, d, e ]
  //     [ 0, 0, f ]

  extern char THETA_DIR[];  // This is the global directory for the theta function files
  
  // Make the filename "ternary_bool_theta__a_b_c_d_e_f__precision.bindata"
  char filename[200];  
  sprintf(filename, "%sternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  THETA_DIR, form[0], form[1], form[2], form[3], form[4], form[5], precision); 
  

  // Compute the length of the array of long:  n = ceil(precision / 32)
  unsigned long n;
  n = precision >> 5;
  if ((precision % 32) != 0)
    n = n + 1;


  /*
  // Make an array to hold the data  -- SUPERFLUOUS since we assume this is already declared!
  unsigned long theta_read[n];
  theta_bin = theta_read;
  */
    

  // -------------------------------------------------------------------

  // Open the file for reading and check it opened correctly
  ifstream filein;
  filein.open(filename, ios::in | ios::binary);

  if (! filein.is_open())
    { cout << "Error opening input file " << filename; exit (1); }

  // Read the array
  filein.read( (char*) theta_read, sizeof(long) * n);

  // Close the file
  filein.close();

}







///////////////////////////////////////////////////////////////////////
// This reads in a (Magma formatted) q-series and translates it into //
// a bitwise boolean theta function formatted as an array of long.   //
// [ Series is assumed to be increasing exponents, formatted like ]  //
// [           2 + 5*q + x^2 + ... + O(q^100)                     ]  //
// Note: We don't allow the coefficient zero.
///////////////////////////////////////////////////////////////////////

//! \note There is a similar routine for reading in a magma theta function, 
//        but this does it incrementally...

void EncodeTernaryThetaBinary(const char* seriesfilename, long form[6]) {

  extern char THETA_DIR[];  // This is the global directory for the theta function files

  
  // Open the seriesfile for reading
  ifstream seriesfile;
  char inputfilename[200];
  sprintf(inputfilename, "%s%s", THETA_DIR, seriesfilename); 
  seriesfile.open(inputfilename, ios::in);
  if (! seriesfile.is_open())
    { cout << "ReadSeries Error: Error opening file"; exit (1); }


  // Open a dummy file for writing
  ofstream outfile;
  char dummyfilename[200];
  cout << "Here 2.1" << endl;  
  sprintf(dummyfilename, "%sdummy_binary_file____%s", THETA_DIR, seriesfilename); 
  cout << "Here 2.2" << endl;
  cout << " dummyfilename = " << dummyfilename << endl;
  outfile.open(dummyfilename, ios::out | ios::trunc | ios::binary);

  if (! outfile.is_open())
    { cout << "Error opening output file " << dummyfilename; exit (1); }


  // Make temporary long to store the current bits
  unsigned long current[1] = { 0 };  
  mpz_class begin;  // This tells the exponent of the first bit.  (Note: sizeof(long) = 4)
  begin = 0;
  
  
  // Initialization
  char c;
  mpz_class num;
  mpz_class pow = -1;
  bool DoneFlag = false;
  
  // Check if we're finished
  while (DoneFlag == false) {
    
    // Check if we're still in the current long
    do {
      
      // Check to see if we're done "O(q^...)"
      seriesfile >> c;
      if (c == 'O') {
	DoneFlag = true;
      }
      
      // Otherwise, try to read the exponent
      else {
	seriesfile.putback(c);
	
	// 1. Check to see of there is no leading coefficient (which means the coefficient is 1)
	if (c == 'q') {
	  num = 1;
	  c = '*';
	}
	else {
	  // Get the coefficient
	  seriesfile >> num;         	
	  
	  // Eat the "+" or "*"
	  seriesfile >> c;
	}
	
	// If there's a '+', then we have a constant term
	if (c == '+') 
	  pow = 0;
	else {
	  
	  // 2. Check if there is a linear term
	  // Eat the "q +" or "q^"
	  seriesfile >> c;
	  seriesfile >> c;
	  
	  // If there's a '+', then we have a linear term
	  if (c == '+') 
	    pow = 1;
	  else {

	    // 3. Otherwise we have a higher-order term
	    // See where to put the term if it has an exponent
	    seriesfile >> pow;
	    
	    // Eat the final "+"
	    seriesfile >> c;
	  }
	}      
	
	// Set the pow bit in the long "current" if it's not too big
	if (pow < (begin + 32)) {
	  /*
	    cout << " pow =  " << pow << endl;
	    cout << " trying to set the bit " << mpz_class(pow - begin).get_ui() << endl;
	    cout << " using the value " << (unsigned long) (1 << mpz_class(pow - begin).get_ui()) << endl;
	    cout << "   current = " << current[0] << endl;
	  */
	  current[0] = current[0] | (1 << mpz_class(pow - begin).get_ui());
	  /*
	    cout << "   current = " << current[0] << endl;
	  */
	}      
	
      }  // End of the  "if (c != 'O')"  clause
      
    } while ((DoneFlag == false) && (pow < begin + 8 * sizeof(long)));
    
    
    // -->  We're done or we're left the current long
    
    // Write the current long, and keep writing zeros until we get back in range
    //    /*
    if (DoneFlag == true) {
      cout << "/n pow = " << pow << endl;
      cout << " begin = " << begin << endl;
      cout << " current long = " << current[0] << endl;
      cout << " sizeof(long) = " << sizeof(long) << endl;
    }
      //    */
    outfile.write( (char*) current, sizeof(long));
    begin = begin + 32;
    current[0] = 0;
    while (begin + 32 <= pow) {
      outfile.write( (char*) current, sizeof(long)); 
      begin = begin + 32;
      cout << " Writing zeros since we're not in the right 4 bytes" << endl; 
      cout << "   begin = " << begin << "   pow = " << pow << endl; 
    }
    
    // Now set the appropriate new bit
    current[0] = current[0] | (1 << (mpz_class(pow - begin).get_ui()));
    
  }

  
  // Find the current precision
  seriesfile >> c;    // Eat the '('
  cout << " just ate: " << c << endl;
  seriesfile >> c;    // Eat the 'q'
  cout << " just ate: " << c << endl;
  seriesfile >> c;    // Eat the '^'
  cout << " just ate: " << c << endl;
  seriesfile >> num;    // Find the precision (num - 1)
  num = num - 1; 
  
  cout << " The precision seems to be " << num << endl;


  // Close all of the files
  seriesfile.close();
  outfile.close();


  // Make a new filename
  char newfilename[200];  
  if (num.fits_ulong_p() == true) {
    sprintf(newfilename, "%sternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	    THETA_DIR, form[0], form[1], form[2], form[3], form[4], form[5], num.get_ui());   
  }
  else {
    cout << "Error in EncodeTernaryThetaBinary: The precision " << num 
	 << " doesn't fit in a long! =(" << endl;
    exit(1);
  }


  // Rename the dummy file 
  int result = rename(dummyfilename, newfilename);
  if (result != 0 ) {
    perror( "\n Error in EncodeTernaryThetaBinary: Renaming the dummy file failed... =( \n\n" );
    exit(1);
  }

  
}






