
PowerSeries<mpz_class> Theta_PARI_1(const Matrix_mpz & Q, const mpz_class & precision) {

  // Sanity check that the precision fits inside an unsigned long
  assert(precision.fits_ulong_p() == true);


  PowerSeries<mpz_class> Theta_Series(precision.get_ui());

  // Note: We assume that pari_init(...) was run at the beginning of the program!


  // Compute the theta function 
  GEN T = qfrep0(Matrix_mpz__to__PARI(Q), stoi(2*precision.get_ui()), 0);

  // Write to thte theta series
  Theta_Series[0] = 1;
  for(unsigned long i = 1; i<=precision.get_ui(); i++)
    Theta_Series[i] =  2 * mpz_class(itos(compo(T,2*i)));  // This assumes we never overflow the long!
 

  return Theta_Series;

}


boolean_theta Theta_PARI_2(const Matrix_mpz & Q, const mpz_class & precision) {


  // Sanity check that the precision fits inside an unsigned long -- so the computation isn't *way* too long!
  assert(precision.fits_ulong_p() == true);


  boolean_theta Theta_Series(Q, precision);

  // Note: We assume that pari_init(...) was run at the beginning of the program!

  /*
  // DIAGNOSTIC
  output(Matrix_mpz__to__PARI(Q));
  output(stoi(precision));
  */

  // Compute the theta function 
  GEN T = qfrep0(Matrix_mpz__to__PARI(Q), stoi(2*precision.get_ui()), 0);

  // Write to the theta series
  Theta_Series.set_value(0);
  for(mpz_class i = 1; i<=precision; i++)
    if (gcmp0(compo(T,2*i.get_ui())) == false)
      Theta_Series.set_value(i);
 

  return Theta_Series;

}


// Compute the theta function by importing (and simplifying) the original PARI code!
PowerSeries<mpz_class> Theta_PARI_1_new(const Matrix_mpz & Q, const mpz_class & precision) {

  // Sanity check that the precision fits inside an unsigned long
  assert(precision.fits_ulong_p() == true);


  // Make the theta series, and set zero.
  PowerSeries<mpz_class> Theta_Series(precision.get_ui());
  Theta_Series[0] = 1;
  
  // Note: We assume that pari_init(...) was run at the beginning of the program!
  extern pari_sp avma;

  // PARI Stack placeholder
  pari_sp av0 = avma;
  
  
  // Make an array of double for the Cholesky decompostion of our matrix
  GEN A = Matrix_mpz__to__PARI(Q);
  
  // Compute the LLL Basis for A
  GEN U = lllgramint(A);
  /*
    // DIAGNSOTIC:
    cout << " lg(U) = " << lg(U) << endl;
    cout << " Q.NumRows() = " << Q.NumRows() << endl;
  */
  if (lg(U) != Q.NumRows()+1) 
    err(talker, "Dude, your matrix is not positive definite! =(");
  
  // Put A in LLL reduced form
  //A = U^T * A * U
  A = qf_base_change(A,U,1);
  
  // Compute the Cholesky decomposition of A
  GEN R = sqred(A);




  // --------------------------------------------


  // Delare some local variables
  long n = Q.NumRows(), i, j, k;
  double p, BOUND, eps = 0.000001;
  BOUND =  2*((double) precision.get_d()) + eps; 


  // Convert R into doubles...
  vector<long> x(n+1);      // Note: This is ok as long since it maxes out at 2 * 10^9, which gives precision at least 4 * 10^18! =)
  vector<double> v(n+1), y(n+1), z(n+1);
  vector< vector <double> > q(n+1);
  for(long i=1; i<=n; i++)
    q[i].resize(n+1);


  /*
  // DIAGNSOTIC:
  cout << "Starting the conversion..." << endl;
  output(R);
  cout << " n = " << n << endl;
  */

  // Set v and q
  for(long i=1; i<=n; i++)
    v[i] = gtodouble(gcoeff(R,i,i));
  for(long i=1; i<=n; i++)
    for(long j=i+1; j<=n; j++) {      
      q[i][j] = gtodouble(gcoeff(R,i,j));
      //cout << "Q[" << i << "][" << j << "] = " << q[i][j] << endl;
    }


  /*
  // DIAGNSOTIC:
  cout << " v = " << v << endl;
  cout << " q = " << q << endl;
  */



  // Reset the PARI stack
  avma = av0;



  /*
  // DIAGNOSTIC
  cout << " Precision = " << precision << endl;  // WARNING: This is a unsigned long, so its biggest value is about 4.2 * 10^9 !
  cout << " BOUND = " << BOUND << endl;
  */



  
  // Set some bounds before starting the main loop
  k = n; y[n] = z[n] = 0;
  x[n] = (long) ceil(sqrt(BOUND/v[n]));

  /*
  // DIAGNOSTIC
  cout << " BOUND/v[n] = " << BOUND/v[n] << endl;
  cout << " sqrt(BOUND/v[n]) = " << sqrt(BOUND/v[n]) << endl;
  cout << " ceil(sqrt(BOUND/v[n])) = " << ceil(sqrt(BOUND/v[n])) << endl;
  cout << " (long) ceil(sqrt(BOUND/v[n])) = " << ((long) ceil(sqrt(BOUND/v[n]))) << endl;
  */

  

  // Compute the theta function
  for(;;x[1]--)
  {
    do
    {
      /*
      cout << " x = " << x << endl;
      */

      // Set more components of x, and the auxilliary y and z.
      if (k>1)
      {
        long l = k-1;

	z[l] = 0;
	for (j=k; j<=n; j++) z[l] += q[l][j]*x[j];

	p = (double) x[k] + z[k];
	y[l] = y[k] + p*p*v[k];

	x[l] = (long) floor(sqrt((BOUND-y[l])/v[l])-z[l]);
        k = l;
      }

      // Check that so far the partial Q(x) is small enough
      for(;;)
      {
	p = (double) x[k] + z[k];
	if (y[k] + p*p*v[k] <= BOUND) break;
	k++; x[k]--;
      }
    }
    while (k > 1);


    // Here we have a valid vector x, with Q(x) <= precision 
    if (! x[1] && y[1]<=eps) break;


    // Compute Q(x)
    p = (double) x[1] + z[1]; 
    p = y[1] + p*p*v[1];   

    // Set the entry in out Theta_Series
    ulong norm = (ulong)(p/2 + 0.5);        // This does the appropriate rounding! =)
    Theta_Series[norm] +=  2;

    // NICE SANITY CHECK FOR LARGE NUMBERS: We could test each p to make sure it's even! =)

  }


  return Theta_Series;

}


// Compute the theta function by importing (and simplifying) the original PARI code!
boolean_theta Theta_PARI_2_new(const Matrix_mpz & Q, const mpz_class & precision) {

  // Make the theta series, and set zero.
  boolean_theta Theta_Series(Q, precision);    // WARNING: This assumes the precision of half of the maximum allowed!  
  Theta_Series.set_value(0);
  
  // Note: We assume that pari_init(...) was run at the beginning of the program!
  extern pari_sp avma;

  // PARI Stack placeholder
  pari_sp av0 = avma;
  
  
  // Make an array of double for the Cholesky decompostion of our matrix
  GEN A = Matrix_mpz__to__PARI(Q);
  
  // Compute the LLL Basis for A
  GEN U = lllgramint(A);
  /*
    // DIAGNSOTIC:
    cout << " lg(U) = " << lg(U) << endl;
    cout << " Q.NumRows() = " << Q.NumRows() << endl;
  */
  if (lg(U) != Q.NumRows()+1) 
    err(talker, "Dude, your matrix is not positive definite! =(");
  
  // Put A in LLL reduced form
  //A = U^T * A * U
  A = qf_base_change(A,U,1);
  
  // Compute the Cholesky decomposition of A
  GEN R = sqred(A);




  // --------------------------------------------


  // Delare some local variables
  long n = Q.NumRows(), i, j, k;
  double p, BOUND, eps = 0.000001;
  BOUND =  2*((double) precision.get_d()) + eps; 


  // Convert R into doubles...
  vector<long> x(n+1);      // Note: This is ok as long since it maxes out at 2 * 10^9, which gives precision at least 4 * 10^18! =)
  vector<double> v(n+1), y(n+1), z(n+1);
  vector< vector <double> > q(n+1);
  for(long i=1; i<=n; i++)
    q[i].resize(n+1);


  /*
  // DIAGNSOTIC:
  cout << "Starting the conversion..." << endl;
  output(R);
  cout << " n = " << n << endl;
  */

  // Set v and q
  for(long i=1; i<=n; i++)
    v[i] = gtodouble(gcoeff(R,i,i));
  for(long i=1; i<=n; i++)
    for(long j=i+1; j<=n; j++) {      
      q[i][j] = gtodouble(gcoeff(R,i,j));
      //cout << "Q[" << i << "][" << j << "] = " << q[i][j] << endl;
    }


  /*
  // DIAGNSOTIC:
  cout << " v = " << v << endl;
  cout << " q = " << q << endl;
  */



  // Reset the PARI stack
  avma = av0;



  /*
  // DIAGNOSTIC
  cout << " Precision = " << precision << endl;  // WARNING: This is a unsigned long, so its biggest value is about 4.2 * 10^9 !
  cout << " BOUND = " << BOUND << endl;
  */



  
  // Set some bounds before starting the main loop
  k = n; y[n] = z[n] = 0;
  x[n] = (long) ceil(sqrt(BOUND/v[n]));

  /*
  // DIAGNOSTIC
  cout << " BOUND/v[n] = " << BOUND/v[n] << endl;
  cout << " sqrt(BOUND/v[n]) = " << sqrt(BOUND/v[n]) << endl;
  cout << " ceil(sqrt(BOUND/v[n])) = " << ceil(sqrt(BOUND/v[n])) << endl;
  cout << " (long) ceil(sqrt(BOUND/v[n])) = " << ((long) ceil(sqrt(BOUND/v[n]))) << endl;
  */

  

  // Compute the theta function
  for(;;x[1]--)
  {
    do
    {
      /*
      cout << " x = " << x << endl;
      */

      // Set more components of x, and the auxilliary y and z.
      if (k>1)
      {
        long l = k-1;

	z[l] = 0;
	for (j=k; j<=n; j++) z[l] += q[l][j]*x[j];

	p = (double) x[k] + z[k];
	y[l] = y[k] + p*p*v[k];

	x[l] = (long) floor(sqrt((BOUND-y[l])/v[l])-z[l]);
        k = l;
      }

      // Check that so far the partial Q(x) is small enough
      for(;;)
      {
	p = (double) x[k] + z[k];
	if (y[k] + p*p*v[k] <= BOUND) break;
	k++; x[k]--;
      }
    }
    while (k > 1);


    // Here we have a valid vector x, with Q(x) <= precision 
    if (! x[1] && y[1]<=eps) break;


    // Compute Q(x)
    p = (double) x[1] + z[1]; 
    p = y[1] + p*p*v[1];   

    // Set the entry in out Theta_Series
    mpz_class norm = mpz_class(p/2 + 0.5);        // This does the appropriate rounding! =)
    Theta_Series.set_value(norm);

    // NICE SANITY CHECK FOR LARGE NUMBERS: We could test each p to make sure it's even! =)

  }


  return Theta_Series;

}




// Compute the theta function by importing (and simplifying) the original PARI code!
boolean_theta Theta_PARI_2_new_Approximate(const Matrix_mpz & Q, const mpz_class & precision, const float & approx_size) {

  // Make the theta series, and set zero.
  boolean_theta Theta_Series(Q, precision);    // WARNING: This assumes the precision of half of the maximum allowed!  
  Theta_Series.set_value(0);
  
  // Note: We assume that pari_init(...) was run at the beginning of the program!
  extern pari_sp avma;

  // PARI Stack placeholder
  pari_sp av0 = avma;
  
  
  // Make an array of double for the Cholesky decompostion of our matrix
  GEN A = Matrix_mpz__to__PARI(Q);
  
  // Compute the LLL Basis for A
  GEN U = lllgramint(A);
  /*
    // DIAGNSOTIC:
    cout << " lg(U) = " << lg(U) << endl;
    cout << " Q.NumRows() = " << Q.NumRows() << endl;
  */
  if (lg(U) != Q.NumRows()+1) 
    err(talker, "Dude, your matrix is not positive definite! =(");
  
  // Put A in LLL reduced form
  //A = U^T * A * U
  A = qf_base_change(A,U,1);
  cout << " The LLL reduction of this matrix is given by:" << endl;
  outmat(A);
  
  // Compute the Cholesky decomposition of A
  GEN R = sqred(A);




  // --------------------------------------------


  // Delare some local variables
  long n = Q.NumRows(), i, j, k;
  double p, BOUND, eps = 0.000001;
  BOUND =  2 * precision.get_d() + eps; 


  // Convert R into doubles...
  vector<long> x(n+1);      // Note: This is ok as long since it maxes out at 2 * 10^9, which gives precision at least 4 * 10^18! =)
  vector<double> v(n+1), y(n+1), z(n+1);
  vector< vector <double> > q(n+1);
  for(long i=1; i<=n; i++)
    q[i].resize(n+1);


  /*
  // DIAGNSOTIC:
  cout << "Starting the conversion..." << endl;
  output(R);
  cout << " n = " << n << endl;
  */

  // Set v and q
  for(long i=1; i<=n; i++)
    v[i] = gtodouble(gcoeff(R,i,i));
  for(long i=1; i<=n; i++)
    for(long j=i+1; j<=n; j++) {      
      q[i][j] = gtodouble(gcoeff(R,i,j));
      //cout << "Q[" << i << "][" << j << "] = " << q[i][j] << endl;
    }


  /*
  // DIAGNSOTIC:
  cout << " v = " << v << endl;
  cout << " q = " << q << endl;
  */



  // Reset the PARI stack
  avma = av0;



  /*
  // DIAGNOSTIC
  cout << " Precision = " << precision << endl;  // WARNING: This is a unsigned long, so its biggest value is about 4.2 * 10^9 !
  cout << " BOUND = " << BOUND << endl;
  */



  // Check for a precision overflow (10^12)
  mpz_class THOUSAND(1000);
  mpz_class PRECISON_OVERFLOW = THOUSAND * THOUSAND * THOUSAND * THOUSAND;
  assert( precision <= PRECISON_OVERFLOW);    // THIS NEEDS TO BE MODIFIED TO DECIDE THE CORRECT CUTOFF!! =)



  
  // Set some bounds before starting the main loop
  float APPROX_SIZE = approx_size;                              // THIS DETERMINES THE AMOUNT WE TRUNCATE IN THE APPROXIMATE ROUTINE! (Default = 500.)
  k = n; y[n] = z[n] = 0;
  x[n] = (long) ceil(sqrt(BOUND/v[n]));

  //   /*
  // DIAGNOSTIC
  cout << " BOUND/v[n] = " << BOUND/v[n] << endl;
  cout << " sqrt(BOUND/v[n]) = " << sqrt(BOUND/v[n]) << endl;
  cout << " ceil(sqrt(BOUND/v[n])) = " << ceil(sqrt(BOUND/v[n])) << endl;
  cout << " (long) ceil(sqrt(BOUND/v[n])) = " << ((long) ceil(sqrt(BOUND/v[n]))) << endl;
  //   */

  

  // Compute the theta function
  for(;;x[1]--)
  {
    // DIAGNOSTIC
    cout << " Using x[1] = " << x[1] << endl;
    cout << " Using x = " << x << endl;
    PrintTime();
    
    do
    {
      /*
      cout << " x = " << x << endl;
      */

      // Set more components of x, and the auxilliary y and z.
      if (k>1)
      {
        long l = k-1;

	z[l] = 0;
	for (j=k; j<=n; j++) z[l] += q[l][j]*x[j];

	p = (double) x[k] + z[k];
	y[l] = y[k] + p*p*v[k];

	x[l] = (long) min( floor(sqrt((BOUND-y[l])/v[l])-z[l]), APPROX_SIZE-z[l] );     // THIS MAKES IT APPROXIMATE FROM ABOVE!!!
        k = l;
      }

      // Check that so far the partial Q(x) is small enough
      for(;;)
      {
	p = (double) x[k] + z[k];
	if ( (k == n) || ( (k != n) && (p >= -APPROX_SIZE) ) )              // THIS MAKES IT APPROXIMATE FROM BELOW!!!
	  if (y[k] + p*p*v[k] <= BOUND) break;
	k++; x[k]--;
      }
      
    }
    while (k > 1);
    

    // Here we have a valid vector x, with Q(x) <= precision 
    if (! x[1] && y[1]<=eps) break;


    // Compute Q(x)
    p = (double) x[1] + z[1]; 
    p = y[1] + p*p*v[1];   


    // SANITY (parity) CHECK 
    // Check to see if the roundoff error is enough to make an odd number.
    // (This should detect roughly half of the roundoff errors!)
    mpz_class even_testnum = mpz_class(p + 0.5);
    assert(even_testnum % 2 == 0);


    // Set the entry in our Theta_Series
    //    cout << " About to set the value for p = " << p << endl; 
    mpz_class norm = mpz_class(p/2 + 0.5);        // This does the appropriate rounding! =)
    //    cout << " Converted p/2 into the number norm = " << norm << endl;
    Theta_Series.set_value(norm);
    /*
    cout << " Finished setting the index for ThetaSeries[" << norm << "]" << endl;
    cout << endl << endl;
    */

    // NICE SANITY CHECK FOR LARGE NUMBERS: We could test each p to make sure it's even! =)

  }


  return Theta_Series;

}






// Compute the theta function by importing (and simplifying) the original PARI code!
boolean_theta Theta_PARI_3_new_Approximate_Ternary(const Matrix_mpz & Q, const mpz_class & precision, const float & approx_size) {

  // Make the theta series, and set zero.
  boolean_theta Theta_Series(Q, precision);    // WARNING: This assumes the precision of half of the maximum allowed!  
  Theta_Series.set_value(0);
  
  // Note: We assume that pari_init(...) was run at the beginning of the program!
  extern pari_sp avma;

  // PARI Stack placeholder
  pari_sp av0 = avma;
  
  
  // Make an array of double for the Cholesky decompostion of our matrix
  GEN A = Matrix_mpz__to__PARI(Q);
  
  // Compute the LLL Basis for A
  GEN U = lllgramint(A);
  /*
    // DIAGNSOTIC:
    cout << " lg(U) = " << lg(U) << endl;
    cout << " Q.NumRows() = " << Q.NumRows() << endl;
  */
  if (lg(U) != Q.NumRows()+1) 
    err(talker, "Dude, your matrix is not positive definite! =(");
  
  // Put A in LLL reduced form
  //A = U^T * A * U
  A = qf_base_change(A,U,1);
  cout << " The LLL reduction of this matrix is given by:" << endl;
  outmat(A);
  
  // Compute the Cholesky decomposition of A
  GEN R = sqred(A);




  // --------------------------------------------


  // Delare some local variables
  long n = Q.NumRows(), i, j, k;
  double p, BOUND, eps = 0.000001;
  double Val;
  BOUND =  2 * precision.get_d() + eps; 


  // Convert R into doubles...
  vector<long> x(n+1);      // Note: This is ok as long since it maxes out at 2 * 10^9, which gives precision at least 4 * 10^18! =)
  vector<double> v(n+1), y(n+1), z(n+1);
  vector< vector <double> > q(n+1);
  for(long i=1; i<=n; i++)
    q[i].resize(n+1);


  /*
  // DIAGNSOTIC:
  cout << "Starting the conversion..." << endl;
  output(R);
  cout << " n = " << n << endl;
  */

  // Set v and q
  for(long i=1; i<=n; i++)
    v[i] = gtodouble(gcoeff(R,i,i));
  for(long i=1; i<=n; i++)
    for(long j=i+1; j<=n; j++) {      
      q[i][j] = gtodouble(gcoeff(R,i,j));
      //cout << "Q[" << i << "][" << j << "] = " << q[i][j] << endl;
    }


  /*
  // DIAGNSOTIC:
  cout << " v = " << v << endl;
  cout << " q = " << q << endl;
  */



  // Reset the PARI stack
  avma = av0;



  /*
  // DIAGNOSTIC
  cout << " Precision = " << precision << endl;  // WARNING: This is a unsigned long, so its biggest value is about 4.2 * 10^9 !
  cout << " BOUND = " << BOUND << endl;
  */



  // Check for a precision overflow (10^12)
  mpz_class THOUSAND(1000);
  mpz_class PRECISON_OVERFLOW = THOUSAND * THOUSAND * THOUSAND * THOUSAND;
  assert( precision <= PRECISON_OVERFLOW);    // THIS NEEDS TO BE MODIFIED TO DECIDE THE CORRECT CUTOFF!! =)



  
  // Set some bounds before starting the main loop
  float APPROX_SIZE = approx_size;   // THIS DETERMINES THE AMOUNT WE TRUNCATE IN THE APPROXIMATE ROUTINE! (Default = 500.)
  //  k = n; y[n] = z[n] = 0;
  //x[n] = (long) ceil(sqrt(BOUND/v[n]));

  //   /*
  // DIAGNOSTIC
  cout << " BOUND/v[n] = " << BOUND/v[n] << endl;
  cout << " sqrt(BOUND/v[n]) = " << sqrt(BOUND/v[n]) << endl;
  cout << " ceil(sqrt(BOUND/v[n])) = " << ceil(sqrt(BOUND/v[n])) << endl;
  cout << " (long) ceil(sqrt(BOUND/v[n])) = " << ((long) ceil(sqrt(BOUND/v[n]))) << endl;
  //   */


  // Sanity Check -- make sure it's a ternary form!
  assert(n==3);

  // Loop through all vectors in the cylinder to make the approximate theta function! =)
  // -----------------------------------------------------------------------------------
  for(x[3] = (long) ceil(sqrt(BOUND/v[3])); x[3] >= 0; x[3]--) {
    
    // Make the cumulative value y[2] from x[3]
    p = (double) x[3];
    y[2] = p*p * v[3];


    for(x[2] = 0; x[2] <= APPROX_SIZE; x[2]++) {

      // Make the cumulative value y[1] from x and y[2]
      p = ((double) x[2])  +  q[2][3] * ((double) x[3]);
      y[1] = p*p * v[2];


      for(x[1] = 0; x[1] <= APPROX_SIZE; x[1]++) {

	// Make the cumulative value Val from x and y[1]
	p = ((double) x[1])  +  q[1][2] * ((double) x[2])   +  q[1][3] * ((double) x[3]);
	Val = p*p * v[1];

	/*
	// DIAGNOSTIC:
	cout << " Using x = " << x << endl;
	cout << " x[1] = " << x[1] << endl;
	cout << " x[2] = " << x[2] << endl;
	cout << " x[3] = " << x[3] << endl;
	cout << " This has value " << Val << endl;
	*/
	
	// SANITY (parity) CHECK 
	// Check to see if the roundoff error is enough to make an odd number.
	// (This should detect roughly half of the roundoff errors!)
	mpz_class even_testnum = mpz_class(Val + 0.5);
	if (even_testnum % 2 != 0) {
	  cout << " Using x = " << x << endl;
	  cout << " x[1] = " << x[1] << endl;
	  cout << " x[2] = " << x[2] << endl;
	  cout << " x[3] = " << x[3] << endl;
	  cout << " This has value " << Val << endl;
	  cout << " even_testnum = " << even_testnum << endl;
	}
	assert(even_testnum % 2 == 0);
	
	
	// Set the entry in our Theta_Series
	//    cout << " About to set the value for p = " << p << endl; 
	mpz_class norm = mpz_class(Val/2 + 0.5);        // This does the appropriate rounding! =)
	//    cout << " Converted p/2 into the number norm = " << norm << endl;
	Theta_Series.set_value(norm);
	/*
	  cout << " Finished setting the index for ThetaSeries[" << norm << "]" << endl;
	  cout << endl << endl;
	*/
	
	// NICE SANITY CHECK FOR LARGE NUMBERS: We could test each p to make sure it's even! =)

      }
	
    }
    
  }


  return Theta_Series;

}







// Compute the theta function by importing (and simplifying) the original PARI code!
boolean_theta Theta_PARI_3_new_Approximate_Ternary_mpz(const Matrix_mpz & Q, const mpz_class & precision, const float & approx_size) {

  // Make the theta series, and set zero.
  boolean_theta Theta_Series(Q, precision);    // WARNING: This assumes the precision of half of the maximum allowed!  
  Theta_Series.set_value(0);
  
  // Note: We assume that pari_init(...) was run at the beginning of the program!
  extern pari_sp avma;

  // PARI Stack placeholder
  pari_sp av0 = avma;
  
  
  // Make an array of double for the Cholesky decompostion of our matrix
  GEN A = Matrix_mpz__to__PARI(Q);
  
  // Compute the LLL Basis for A
  GEN U = lllgramint(A);
  /*
    // DIAGNSOTIC:
    cout << " lg(U) = " << lg(U) << endl;
    cout << " Q.NumRows() = " << Q.NumRows() << endl;
  */
  if (lg(U) != Q.NumRows()+1) 
    err(talker, "Dude, your matrix is not positive definite! =(");
  
  // Put A in LLL reduced form
  //A = U^T * A * U
  A = qf_base_change(A,U,1);
  cout << " The LLL reduction of this matrix is given by:" << endl;
  outmat(A);


  // Recreate the matrix
  mpz_class q11(gtolong(gcoeff(A,1,1)));
  mpz_class q22(gtolong(gcoeff(A,2,2)));
  mpz_class q33(gtolong(gcoeff(A,3,3)));

  mpz_class q12(2 * gtolong(gcoeff(A,1,2)));
  mpz_class q13(2 * gtolong(gcoeff(A,1,3)));
  mpz_class q23(2 * gtolong(gcoeff(A,2,3)));


  // DIAGNOSTIC
  cout << " We now have the quadratic form " 
       << q11 << "x^2 + "
       << q22 << "y^2 + "
       << q33 << "z^2 + "
       << q12 << "xy + "
       << q13 << "xz + "
       << q23 << "yz." << endl;

    

  // Delare some local variables
  long n = Q.NumRows(), i, j, k;
  double BOUND =  2 * precision.get_d() + 0.0001; 
  vector<mpz_class> x(4); 
  mpz_class Val, y2, y3;


  // Compute the maximum x[3]-value we will need
  mpz_class ternary_disc(4*q11*q22*q33 - (q11*q23 + q22*q13 + q33*q12) + q12*q23*q13);
  mpz_class x3_disc(4*q11*q22 - q12*q12);
  long x3_max = (long) ceil(sqrt(BOUND * x3_disc.get_d() / ternary_disc.get_d()));
  long old_x3_max = (long) ceil(sqrt(BOUND/q33.get_d()));
				
				
  // DIAGNOSTIC
  cout << " The previous (very naive) w-bound was: " << old_x3_max << endl;
  cout << " The new (correct) w-bound is: " << x3_max << endl;


  // Reset the PARI stack
  avma = av0;



  /*
  // DIAGNOSTIC
  cout << " Precision = " << precision << endl;  // WARNING: This is a unsigned long, so its biggest value is about 4.2 * 10^9 !
  cout << " BOUND = " << BOUND << endl;
  */



  // Check for a precision overflow (10^12)
  mpz_class THOUSAND(1000);
  mpz_class PRECISON_OVERFLOW = THOUSAND * THOUSAND * THOUSAND * THOUSAND;
  assert( precision <= PRECISON_OVERFLOW);    // THIS NEEDS TO BE MODIFIED TO DECIDE THE CORRECT CUTOFF!! =)



  
  // Set some bounds before starting the main loop
  float APPROX_SIZE = approx_size;   // THIS DETERMINES THE AMOUNT WE TRUNCATE IN THE APPROXIMATE ROUTINE! (Default = 500.)
  //  k = n; y[n] = z[n] = 0;
  //x[n] = (long) ceil(sqrt(BOUND/v[n]));

  /*
  // DIAGNOSTIC
  cout << " BOUND/v[n] = " << (BOUND/q33.get_d()) << endl;
  cout << " sqrt(BOUND/v[n]) = " << sqrt(BOUND/q33.get_d()) << endl;
  cout << " ceil(sqrt(BOUND/v[n])) = " << ceil(sqrt(BOUND/q33.get_d())) << endl;
  cout << " (long) ceil(sqrt(BOUND/v[n])) = " << ((long) ceil(sqrt(BOUND/q33.get_d()))) << endl;
  */


  // Sanity Check -- make sure it's a ternary form!
  assert(n==3);

  // Loop through all vectors in the cylinder to make the approximate theta function! =)
  // -----------------------------------------------------------------------------------
  for(x[3] = x3_max; x[3] >= 0; x[3]--) {
    
    // Make the cumulative value y3 from x
    y3 = q33 * x[3] * x[3];


    for(x[2] = 0; x[2] <= APPROX_SIZE; x[2]++) {

      // Make the cumulative value y2 from x and y3
      y2 = q22*x[2]*x[2]  + q23*x[2]*x[3] + y3;


      for(x[1] = 0; x[1] <= APPROX_SIZE; x[1]++) {

	// Make the cumulative value Val from x and y2
	Val = q11*x[1]*x[1]  + q12*x[1]*x[2] + q13*x[1]*x[3] + y2;


	/*
	// DIAGNOSTIC:
	cout << " Using x = " << x << endl;
	cout << " x[1] = " << x[1] << endl;
	cout << " x[2] = " << x[2] << endl;
	cout << " x[3] = " << x[3] << endl;
	cout << " This has value " << Val << endl;
	*/
	
	// SANITY (parity) CHECK 
	// Check to see if the roundoff error is enough to make an odd number.
	// (This should detect roughly half of the roundoff errors!)
	if (Val % 2 != 0) {
	  cout << "WARNING: Parity check failure! =0" << endl;
	  cout << " Using x = " << x << endl;
	  cout << " x[1] = " << x[1] << endl;
	  cout << " x[2] = " << x[2] << endl;
	  cout << " x[3] = " << x[3] << endl;
	  cout << " This has value " << Val << endl;
	}
	assert(Val % 2 == 0);
	
	
	// Set the entry in our Theta_Series

	if (Val <= BOUND) {
	  //    cout << " About to set the value for p = " << p << endl; 
	  mpz_class norm = mpz_class(Val/2);   
	  //    cout << " Converted p/2 into the number norm = " << norm << endl;
	  Theta_Series.set_value(norm);
	  /*
	    cout << " Finished setting the index for ThetaSeries[" << norm << "]" << endl;
	    cout << endl << endl;
	  */
	}
	
	// NICE SANITY CHECK FOR LARGE NUMBERS: We could test each p to make sure it's even! =)

      }
	
    }
    
  }


  return Theta_Series;

}








// ============================================== Slow Theta Function Routines =============================================================



PowerSeries<mpz_class> Theta1(const Matrix_mpz & Q, const mpz_class & precision) {

  // Sanity check that the precision fits inside an unsigned long
  assert(precision.fits_ulong_p() == true);


  PowerSeries<mpz_class> Theta_Series(precision.get_ui());


  
  // Define the value at zero to be -1, to avoid testing it all of the time! =)
  Theta_Series[0] = -1;


  
  /*
  for(long i=0; i < min(100, precision); i++)
    cout << "i = " << i << "    Theta[i] = " << Theta_Series[i] << endl;
  */



  cout << " This will compute the theta series by counting points quickly in a box! =) " << endl; 

  mpz_class desired_precision = precision;


  // Basic Sanity checks
  assert(Q.IsSymmetric() == true);
  assert(Q.NumRows() == 3);


  // For convenience, we define              ****** are these right???  don't we need half of a, b, and c??? ******
  mpz_class a = Q(1,1)/2;
  mpz_class b = Q(2,2)/2;
  mpz_class c = Q(3,3)/2;
  mpz_class d = Q(2,3);
  mpz_class e = Q(1,3);
  mpz_class f = Q(1,2);

  mpz_class disc = 4*a*b*c + d*e*f - a*d*d - b*e*e - c*f*f + d*e*f;

  mpz_class w11 = 4*b*c - d*d;
  mpz_class w22 = 4*a*c - e*e;
  mpz_class w33 = 4*a*b - f*f;


  /*
  // DIAGNOSTIC
  cout << " a = " << a << endl;
  cout << " b = " << b << endl;
  cout << " c = " << c << endl;
  cout << " d = " << d << endl;
  cout << " e = " << e << endl;
  cout << " f = " << f << endl;
  cout << endl;
  cout << " disc = " << disc << endl;
  cout << " w11 = " << w11 << endl;
  cout << " w22 = " << w22 << endl;
  cout << " w33 = " << w33 << endl;
  cout << endl;
  */






  // Set variables for the range of numbers to check:
  // -------------------------------------------------  
  
  //  mpz_class Bunch_Size = 1000000;   // One million
  //  mpz_class Bunch_Size = 100000;   // One Hundred Thousand
  mpz_class Bunch_Size = 10000;   // Ten Thousand
  mpz_class Save_Increment = 100000; // Save the computation after every hundred thousand numbers.

  mpz_class Checking_Min = 0;
  mpz_class Save_Target = Checking_Min + Save_Increment;  

  
  set<mpz_class> eligible_set;


  // REPLACE THIS WITH SOMETHING THAT WORKS ON BUNCHES!!
  mpz_class Checking_Max = desired_precision;

        
    
    // Compute the bounds on vectors to check:
    // ---------------------------------------
    cout << "   Checking possible (locally represented) exceptions from " << Checking_Min << " to " << Checking_Max << " ---- ";

    
    // Approach #1: From Jagy
    mpz_class x_Max, y_Max, z_Max;
    x_Max = ceil(sqrt(abs(Checking_Max.get_d() * w11.get_d() / disc.get_d())));
    y_Max = ceil(sqrt(abs(Checking_Max.get_d() * w22.get_d() / disc.get_d())));
    z_Max = ceil(sqrt(abs(Checking_Max.get_d() * w33.get_d() / disc.get_d())));

    //    /*  
    // DIAGNOSTIC 
    cout << " w11 = " << w11 << "  " << w11.get_d() << endl;
    cout << " w22 = " << w22 << "  " << w22.get_d() << endl;
    cout << " w33 = " << w33 << "  " << w33.get_d() << endl;
    cout << " disc = " << disc << "  " << disc.get_d() << endl;
    //  cout << " precision = " << precision << "  " << precision.get_d() << endl;
    //  cout << " (precision.get_d() * w11.get_d() / disc.get_d()) = " << (precision.get_d() * w11.get_d() / disc.get_d()) << endl;
    cout << " Checking_Max = " << Checking_Max << "  " << Checking_Max.get_d() << endl;
    cout << " (Checking_Max.get_d() * w11.get_d() / disc.get_d()) = " << (Checking_Max.get_d() * w11.get_d() / disc.get_d()) << endl;
    cout << endl;
    cout << " Using the bounds: " << endl;
    cout << " x_Max = " << x_Max << endl;
    cout << " y_Max = " << y_Max << endl;
    cout << " z_Max = " << z_Max << endl;
    cout << endl;  
    // */
  
    
    
    // Cross out the eligible numbers
    // ------------------------------
    
    /*
    // Approach #1: Blindly cross out all numbers in an octant (box)
    mpz_class last_size = eligible_set.size();
    mpz_class diag, D, E, F;
    for(mpz_class x=0; (x <= x_Max) && (eligible_set.empty() == false); x++)
    for(mpz_class y=0; (y <= y_Max) && (eligible_set.empty() == false); y++)
    for(mpz_class z=0; (z <= z_Max) && (eligible_set.empty() == false); z++) {
    
    // Compute the diagonal and off-diagonal terms
    diag = a*x*x + b*y*y + c*z*z;
    D = d*y*z;
    E = e*x*z;
    F = f*x*y;
    
    // Compute and remove the 4 values
    eligible_set.erase(diag + D + E + F);
    eligible_set.erase(diag - D - E + F);   // switch z -> -z
    eligible_set.erase(diag - D + E - F);   // switch y -> -y
    eligible_set.erase(diag + D - E - F);   // switch x -> -x
    
    
    // Check if the size has dropped by 10,000
    if (eligible_set.size() < last_size - 10000) {
    cout << " The current eligible exception set has size " << eligible_set.size() << endl;
    last_size = eligible_set.size();
    }
    
    }	
    */
    
    
    
    // Approach #2: From Jagy's program
    mpz_class TRUANT, MAX_TARGET = Checking_Max;                       // <--- NEED TO SET THESE! =)
    
    mpz_class dy_ex, axx_fxy_byy;
    mpz_class Big_Disc, Small_Disc;
    mpz_class Sqrt_Big_Disc, Sqrt_Small_Disc;
    mpz_class z1, z2, z3, z4;
    
    mpz_class last_size = 0;
    mpz_class diag, D, E, F;
    
    for(mpz_class x=0; (x <= x_Max); x++) 
      
      for(mpz_class y=y_Max; (y >= 0); y--) {
	
	dy_ex = d*y + e*x;                         // Temporary variable
	axx_fxy_byy = a*x*x + f*x*y + b*y*y;       // Temporary variable
	
	// Compute Big_Disc, Small_Disc, and TRUANT
	Big_Disc = (dy_ex) * (dy_ex)  - 4 * c * (axx_fxy_byy - MAX_TARGET);
	Small_Disc = (dy_ex) * (dy_ex)  - 4 * c * (axx_fxy_byy - TRUANT);
	TRUANT = 0;
	
	// Compute the bounds z1 -> z4
	if (Big_Disc >= 0) {
	  Sqrt_Big_Disc = _Jagy_IntSqrt(Big_Disc);         // Temporary Variable
	  
	  z4 = ( -(dy_ex) + Sqrt_Big_Disc ) / (2*c);
	  z1 = ( -(dy_ex) - Sqrt_Big_Disc ) / (2*c);
	  
	  z3 = -(dy_ex);
	  z2 = z3;
	  
	  if (Small_Disc >= 0) {	
	    Sqrt_Small_Disc = _Jagy_IntSqrt(Small_Disc);   // Temporary Variable
	    z3 = ( -(dy_ex) + Sqrt_Small_Disc ) / (2*c);
	    z2 = ( -(dy_ex) - Sqrt_Small_Disc ) / (2*c);
	  }
	  
	}
	
      
	// 2 Loops for z:  z3 --> z4  and  z1 --> z2 
	// ------------------------------------------
	for(mpz_class z=z1; (z <= z4) ; z++) {
	  
	  // Compute the diagonal and off-diagonal terms                           // <--- Can speed this up! =)
	  diag = a*x*x + b*y*y + c*z*z;
	  D = d*y*z;
	  E = e*x*z;
	  F = f*x*y;


	  /*	  
	  // DIAGNOSTIC:
	  mpz_class test_num(25);
	  if ((diag + D + E + F == test_num) || (diag - D - E + F == test_num) || 
	      (diag - D + E - F == test_num) || (diag + D - E - F == test_num)) {
	    cout << endl;
	    cout << " Found the test number " << test_num << endl;
	    cout << "   diag + D + E + F = " << (diag + D + E + F) << endl;
	    cout << "   diag - D - E + F = " << (diag - D - E + F) << endl;
	    cout << "   diag - D + E - F = " << (diag - D + E - F) << endl;
	    cout << "   diag + D - E - F = " << (diag + D - E - F) << endl;
	    cout << endl;	      
	    cout << "   x = " << x << "   y = " << y << "   z = " << z << endl;
	    cout << "   a = " << a << "   b = " << b << "   c = " << c << "   d = " << d << "   e = " << e << "   f = " << f << endl; 
	    cout << "   diag = " << diag << "   D = " << D << "   E = " << E << "   F = " << F << endl;
	    cout << endl;	      	      
	  }
	  */


	  
	  // Compute and remove the 4 values
	  if ((diag + D + E + F) <= precision)
	    Theta_Series[mpz_class(diag + D + E + F).get_ui()] += 2;
	  if ((F != 0) && ((diag - D - E + F) <= precision))
	    Theta_Series[mpz_class(diag - D - E + F).get_ui()] += 2;
	  if ((E != 0) && ((diag - D + E - F) <= precision))
	    Theta_Series[mpz_class(diag - D + E - F).get_ui()] += 2;
	  if ((D != 0) && ((diag + D - E - F) <= precision))
	    Theta_Series[mpz_class(diag + D - E - F).get_ui()] += 2;
	  /*
	  eligible_set.erase(diag + D + E + F);
	  eligible_set.erase(diag - D - E + F);   // switch z -> -z
	  eligible_set.erase(diag - D + E - F);   // switch y -> -y
	  eligible_set.erase(diag + D - E - F);   // switch x -> -x
	  */	  

	}
	
      }
    


  

    /*
  for(long i=0; i < min(100, precision); i++)
    cout << "i = " << i << "    Theta[i] = " << Theta_Series[i] << endl;
    */
  

  return Theta_Series;




}





boolean_theta Theta2(const Matrix_mpz & Q, const mpz_class & precision) {

  // Sanity check that the precision fits inside an unsigned long -- so the computation doesn't take *way* too long!
  assert(precision.fits_ulong_p() == true);


  boolean_theta Theta_Series(Q, precision);


  
  /*
  for(long i=0; i < min(100, precision); i++)
    cout << "i = " << i << "    Theta[i] = " << Theta_Series[i] << endl;
  */



  cout << " This will compute the theta series by counting points quickly in a box! =) " << endl; 

  mpz_class desired_precision = precision;


  // Basic Sanity checks
  assert(Q.IsSymmetric() == true);
  assert(Q.NumRows() == 3);


  // For convenience, we define              ****** are these right???  don't we need half of a, b, and c??? ******
  mpz_class a = Q(1,1)/2;
  mpz_class b = Q(2,2)/2;
  mpz_class c = Q(3,3)/2;
  mpz_class d = Q(2,3);
  mpz_class e = Q(1,3);
  mpz_class f = Q(1,2);

  mpz_class disc = 4*a*b*c + d*e*f - a*d*d - b*e*e - c*f*f + d*e*f;

  mpz_class w11 = 4*b*c - d*d;
  mpz_class w22 = 4*a*c - e*e;
  mpz_class w33 = 4*a*b - f*f;


  /*
  // DIAGNOSTIC
  cout << " a = " << a << endl;
  cout << " b = " << b << endl;
  cout << " c = " << c << endl;
  cout << " d = " << d << endl;
  cout << " e = " << e << endl;
  cout << " f = " << f << endl;
  cout << endl;
  cout << " disc = " << disc << endl;
  cout << " w11 = " << w11 << endl;
  cout << " w22 = " << w22 << endl;
  cout << " w33 = " << w33 << endl;
  cout << endl;
  */






  // Set variables for the range of numbers to check:
  // -------------------------------------------------  
  
  //  mpz_class Bunch_Size = 1000000;   // One million
  //  mpz_class Bunch_Size = 100000;   // One Hundred Thousand
  mpz_class Bunch_Size = 10000;   // Ten Thousand
  mpz_class Save_Increment = 100000; // Save the computation after every hundred thousand numbers.

  mpz_class Checking_Min = 0;
  mpz_class Save_Target = Checking_Min + Save_Increment;  

  
  set<mpz_class> eligible_set;


  // REPLACE THIS WITH SOMETHING THAT WORKS ON BUNCHES!!
  mpz_class Checking_Max = desired_precision;

        
    
    // Compute the bounds on vectors to check:
    // ---------------------------------------
    cout << "   Checking possible (locally represented) exceptions from " << Checking_Min << " to " << Checking_Max << " ---- ";

    
    // Approach #1: From Jagy
    mpz_class x_Max, y_Max, z_Max;
    x_Max = ceil(sqrt(abs(Checking_Max.get_d() * w11.get_d() / disc.get_d())));
    y_Max = ceil(sqrt(abs(Checking_Max.get_d() * w22.get_d() / disc.get_d())));
    z_Max = ceil(sqrt(abs(Checking_Max.get_d() * w33.get_d() / disc.get_d())));

    //    /*  
    // DIAGNOSTIC 
    cout << " w11 = " << w11 << "  " << w11.get_d() << endl;
    cout << " w22 = " << w22 << "  " << w22.get_d() << endl;
    cout << " w33 = " << w33 << "  " << w33.get_d() << endl;
    cout << " disc = " << disc << "  " << disc.get_d() << endl;
    //  cout << " precision = " << precision << "  " << precision.get_d() << endl;
    //  cout << " (precision.get_d() * w11.get_d() / disc.get_d()) = " << (precision.get_d() * w11.get_d() / disc.get_d()) << endl;
    cout << " Checking_Max = " << Checking_Max << "  " << Checking_Max.get_d() << endl;
    cout << " (Checking_Max.get_d() * w11.get_d() / disc.get_d()) = " << (Checking_Max.get_d() * w11.get_d() / disc.get_d()) << endl;
    cout << endl;
    cout << " Using the bounds: " << endl;
    cout << " x_Max = " << x_Max << endl;
    cout << " y_Max = " << y_Max << endl;
    cout << " z_Max = " << z_Max << endl;
    cout << endl;  
    // */
  
    
    
    // Cross out the eligible numbers
    // ------------------------------
    
    /*
    // Approach #1: Blindly cross out all numbers in an octant (box)
    mpz_class last_size = eligible_set.size();
    mpz_class diag, D, E, F;
    for(mpz_class x=0; (x <= x_Max) && (eligible_set.empty() == false); x++)
    for(mpz_class y=0; (y <= y_Max) && (eligible_set.empty() == false); y++)
    for(mpz_class z=0; (z <= z_Max) && (eligible_set.empty() == false); z++) {
    
    // Compute the diagonal and off-diagonal terms
    diag = a*x*x + b*y*y + c*z*z;
    D = d*y*z;
    E = e*x*z;
    F = f*x*y;
    
    // Compute and remove the 4 values
    eligible_set.erase(diag + D + E + F);
    eligible_set.erase(diag - D - E + F);   // switch z -> -z
    eligible_set.erase(diag - D + E - F);   // switch y -> -y
    eligible_set.erase(diag + D - E - F);   // switch x -> -x
    
    
    // Check if the size has dropped by 10,000
    if (eligible_set.size() < last_size - 10000) {
    cout << " The current eligible exception set has size " << eligible_set.size() << endl;
    last_size = eligible_set.size();
    }
    
    }	
    */
    
    
    
    // Approach #2: From Jagy's program
    mpz_class TRUANT, MAX_TARGET = Checking_Max;                       // <--- NEED TO SET THESE! =)
    
    mpz_class dy_ex, axx_fxy_byy;
    mpz_class Big_Disc, Small_Disc;
    mpz_class Sqrt_Big_Disc, Sqrt_Small_Disc;
    mpz_class z1, z2, z3, z4;
    
    mpz_class last_size = 0;
    mpz_class diag, D, E, F;
    
    for(mpz_class x=0; (x <= x_Max); x++) 
      
      for(mpz_class y=y_Max; (y >= 0); y--) {
	
	dy_ex = d*y + e*x;                         // Temporary variable
	axx_fxy_byy = a*x*x + f*x*y + b*y*y;       // Temporary variable
	
	// Compute Big_Disc, Small_Disc, and TRUANT
	Big_Disc = (dy_ex) * (dy_ex)  - 4 * c * (axx_fxy_byy - MAX_TARGET);
	Small_Disc = (dy_ex) * (dy_ex)  - 4 * c * (axx_fxy_byy - TRUANT);
	TRUANT = 0;
	
	// Compute the bounds z1 -> z4
	if (Big_Disc >= 0) {
	  Sqrt_Big_Disc = _Jagy_IntSqrt(Big_Disc);         // Temporary Variable
	  
	  z4 = ( -(dy_ex) + Sqrt_Big_Disc ) / (2*c);
	  z1 = ( -(dy_ex) - Sqrt_Big_Disc ) / (2*c);
	  
	  z3 = -(dy_ex);
	  z2 = z3;
	  
	  if (Small_Disc >= 0) {	
	    Sqrt_Small_Disc = _Jagy_IntSqrt(Small_Disc);   // Temporary Variable
	    z3 = ( -(dy_ex) + Sqrt_Small_Disc ) / (2*c);
	    z2 = ( -(dy_ex) - Sqrt_Small_Disc ) / (2*c);
	  }
	  
	}
	
      
	// 2 Loops for z:  z3 --> z4  and  z1 --> z2 
	// ------------------------------------------
	for(mpz_class z=z1; (z <= z4) ; z++) {
	  
	  // Compute the diagonal and off-diagonal terms                           // <--- Can speed this up! =)
	  diag = a*x*x + b*y*y + c*z*z;
	  D = d*y*z;
	  E = e*x*z;
	  F = f*x*y;


	  /*	  
	  // DIAGNOSTIC:
	  mpz_class test_num(25);
	  if ((diag + D + E + F == test_num) || (diag - D - E + F == test_num) || 
	      (diag - D + E - F == test_num) || (diag + D - E - F == test_num)) {
	    cout << endl;
	    cout << " Found the test number " << test_num << endl;
	    cout << "   diag + D + E + F = " << (diag + D + E + F) << endl;
	    cout << "   diag - D - E + F = " << (diag - D - E + F) << endl;
	    cout << "   diag - D + E - F = " << (diag - D + E - F) << endl;
	    cout << "   diag + D - E - F = " << (diag + D - E - F) << endl;
	    cout << endl;	      
	    cout << "   x = " << x << "   y = " << y << "   z = " << z << endl;
	    cout << "   a = " << a << "   b = " << b << "   c = " << c << "   d = " << d << "   e = " << e << "   f = " << f << endl; 
	    cout << "   diag = " << diag << "   D = " << D << "   E = " << E << "   F = " << F << endl;
	    cout << endl;	      	      
	  }
	  */


	  
	  // Compute and remove the 4 values
	  if ((diag + D + E + F) <= precision)
	    Theta_Series.set_value(mpz_class(diag + D + E + F).get_ui());
	  if ((F != 0) && ((diag - D - E + F) <= precision))
	    Theta_Series.set_value(mpz_class(diag - D - E + F).get_ui());
	  if ((E != 0) && ((diag - D + E - F) <= precision))
	    Theta_Series.set_value(mpz_class(diag - D + E - F).get_ui());
	  if ((D != 0) && ((diag + D - E - F) <= precision))
	    Theta_Series.set_value(mpz_class(diag + D - E - F).get_ui());
	  /*
	  eligible_set.erase(diag + D + E + F);
	  eligible_set.erase(diag - D - E + F);   // switch z -> -z
	  eligible_set.erase(diag - D + E - F);   // switch y -> -y
	  eligible_set.erase(diag + D - E - F);   // switch x -> -x
	  */	  

	}
	
      }
    


  

    /*
  for(long i=0; i < min(100, precision); i++)
    cout << "i = " << i << "    Theta[i] = " << Theta_Series[i] << endl;
    */
  

  return Theta_Series;




}






// Jagy's IntSqrt routine used for his ternary exception finder. =)
// (This returns the largest integer x so that x^2 < n.)                    <=====  CHECK THIS... the end is confusing...
inline 
mpz_class _Jagy_IntSqrt(const mpz_class & n) {
  
  // Check if x=1
  if (n<=0) 
    return 0;
  else if (n<=3)
    return 1;
  
  else {
    mpz_class oldroot = 1;
    mpz_class root = n;
    
    // Do Newton's method to find an approximate root??
    while (abs(root-oldroot) > 5) {
      oldroot = root;
      root = ((n / root) + root) / 2;
    }
    
    // Find the
    while (root * root < n)
      root++;

    while (root * root > n)
      root--;
    
    return root;
  }
  
}




