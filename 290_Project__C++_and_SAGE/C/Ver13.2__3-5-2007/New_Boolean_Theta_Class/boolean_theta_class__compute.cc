

//////////////////////////////////////////////////////////////////////
// Writes the theta function of a ternary form of desired precision //
//////////////////////////////////////////////////////////////////////

void boolean_theta::compute() {
  /*
  // Print the ternary form and its level
  cout << endl;
  cout << " The form is:  [ " 
       << QQ[0] << ", "
       << QQ[1] << ", "
       << QQ[2] << ", "
       << QQ[3] << ", "
       << QQ[4] << ", "
       << QQ[5] << " ] " << endl;
  cout << " The level of the form is " << QF_Ternary_Level(QQ) << endl;
  */


  const long n = 3;


  // Find the (lower-triangular) Cholesky Decomposition (uses indices 1 --> n)
  double* Cholesky;
  CholeskyDecomposition(QQ, Cholesky);

  /*  
  cout << " Computing the theta function " << endl;
  PrintTime();

  cout << endl << "Entering FastBinaryTheta" << endl;
  cout << " Using precision = " << precision() << endl;
  */

  // ERROR CHECKING: Check that C+1 fits in an unsigned long.


  // Zero out the theta function
  for (unsigned long i=0; i < (unsigned long) _length(); i++)
    _theta[i] = 0;


  // Make the (constant) matrix Q from QQ  -- NOTE: We only use indices 1 --> n.
  double Q[n+1][n+1];
  for (long i=1; i<=n; i++) 
    for (long j=1; j<=n; j++) 
      Q[i][j] = 0;   // Clear the matrix
  long counter = 0;
  for (long i=1; i<=n; i++) 
    for (long j=i; j<=n; j++) {
      Q[i][j] = Cholesky[counter];  // Put Cholesky in the upper triangular part
      counter++;
    }


  /*
  // Print Q  
  cout << "Using Q = " << endl;
  for (long i=1; i<=n; i++) {
    for (long j=1; j<=n; j++)
      cout << Q[i][j] << " ";
    cout << endl;
  }
  */


  // 1. Initialize
  long i = n;
  vector<double> T(n+1, 0);  // Note: We use n+1 so we can index the entries as 1 --> n
  vector<double> U(n+1, 0);
  T[i] = (double) precision().get_d();
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
    Q_val_double = precision().get_d() - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    Q_val = (unsigned long) round(Q_val_double);

    //    cout << " Float = " << Q_val_double << "   Long = " << Q_val << "  XX " << endl;
    /*
    cout << " The float value is " << Q_val_double << endl;
    cout << " The associated long value is " << Q_val << endl;
    cout << endl;
    */

    if (Q_val <= precision()) {
      set_value(Q_val);         // Set the appropriate bit 
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


  /*   
  cout << "Leaving FastBinaryTheta" << endl << endl;
  */




  /*
  // DIAGNOSTIC:
   cout << " The last two longs are: " << endl;
  cout << "  i = " << ((precision() >> 5) - 1) << " theta[i] = " << _theta[(precision() >> 5) - 1] << endl;
  cout << "  i = " << ((precision() >> 5) + 0) << " theta[i] = " << _theta[(precision() >> 5) + 0] << endl;
  cout << endl;
  */

}
