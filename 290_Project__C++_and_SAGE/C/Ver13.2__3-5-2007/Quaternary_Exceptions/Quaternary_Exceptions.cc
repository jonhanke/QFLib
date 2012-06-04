

// Initializes the class and computes the exceptions up to a given precision.
QuaternaryExceptions::QuaternaryExceptions(const Matrix_mpz & T, const mpz_class & precision, 
				     const string & quaternary_exceptions_dir) {

  // DIAGNOSTIC
  cout << " Entering QuaternaryExceptions Constructor" << endl;
  cout << " precision = " << precision << endl;
  cout << " quaternary form = " << endl << T << endl;
  

  // Basic Sanity check
  assert( T.NumRows() == T.NumCols() );
  assert( T.NumRows() == 4 );

  // Sanity check for even diagonal (since it should be the matrix of 2*T)
  assert( T(1,1) % 2 == 0 );
  assert( T(2,2) % 2 == 0 );
  assert( T(3,3) % 2 == 0 );
  assert( T(4,4) % 2 == 0 );


  // Set the form
  _T = T;

  // Compute its local conditions
  _local_conditions = LocalConditions(T);


  // DIAGNOSTIC
  cout << " Finished the local conditions, now looking for exceptions! " << endl;

  _QuaternaryExceptionsUpTo(precision, quaternary_exceptions_dir);

}   




// Computes the exceptions from a given set of possibilities
QuaternaryExceptions::QuaternaryExceptions(const Matrix_mpz & T, const set<mpz_class> & possible_set) {
  // DIAGNOSTIC
  cout << " Entering QuaternaryExceptions Constructor" << endl;
  cout << " possible set = " << possible_set << endl;
  cout << " quaternary form = " << endl << T << endl;
  

  // Basic Sanity check
  assert( T.NumRows() == T.NumCols() );
  assert( T.NumRows() == 4 );

  // Sanity check for even diagonal (since it should be the matrix of 2*T)
  assert( T(1,1) % 2 == 0 );
  assert( T(2,2) % 2 == 0 );
  assert( T(3,3) % 2 == 0 );
  assert( T(4,4) % 2 == 0 );


  // Set the form
  _T = T;

  // Compute its local conditions
  _local_conditions = LocalConditions(T);


  // DIAGNOSTIC
  cout << " Finished the local conditions, now looking for exceptions! " << endl;


  // Check each number in the set
  _starting_set = possible_set;
  for(set<mpz_class>::iterator i = possible_set.begin(); i!= possible_set.end(); i++)
  if (_IsQuaternaryException(*i) == true)
    _exceptional_subset.insert(*i);


  // TO CHANGE LATER...
  _exception_set = _exceptional_subset;

}   











// Give the precision x the exceptions are known up to (i.e. all numbers <= x).
mpz_class QuaternaryExceptions::GetPrecision() const {
  return _precision;
}


// Give the set of exceptions (up to _precision).
set<mpz_class> QuaternaryExceptions::GetExceptions() const {
  return _exception_set;
}


// Verify the exceptions (up to the smaller of max and _precision).
bool QuaternaryExceptions::VerifyExceptions(const mpz_class & Max) const {

  // Find the desired number to check up to
  unsigned long Theta_precision;
  if (Max == 0) 
    Theta_precision = _precision.get_ui();            // Check all numbers if no argument is given,
  else 
    Theta_precision = min(Max, _precision).get_ui();  // Otherwise take the minimum of the two.


  // Compute the theta function of _T
  PowerSeries<mpz_class> theta;
  theta = _T.GetMagmaThetaSeries("", Theta_precision);

  // Check that every locally represented number is either represented or in _exception_set
  // (update output every 100,000 numbers checked, and abort on any mistakes)
  for(unsigned long num=0; num <= Theta_precision; num++) {
    if (_local_conditions.IsLocallyRepresented(num) == true) 

      // Check for discrepancies
      if ( ( (_exception_set.find(num) != _exception_set.end()) && (theta[num] > 0) )         // in exception_set but not an exception
	   || ( (_exception_set.find(num) == _exception_set.end()) && (theta[num] == 0) ) ) { // not in exception_set but an exception

	// Print Error message
	cout << " Error in QuaternaryExceptions::VerifyExceptions():  There's a problem with " << num << endl;
	cout << "   " << num << " is in the exception set:  " << (_exception_set.find(num) != _exception_set.end()) << endl;
	cout << "   Theta coefficient at " << num << " = " << theta[num] << " > 0 " << endl;

	// Exit
	assert(0==1);
	return false;
      }	

    // Output when we hit a multiple of 100,000
    if ((num % 100000) == 0)
      cout << "  Ok up to " << num << endl;
    
  }


  // Everything is ok. =)
  return true;

}





// ------------------------------------------------------------------------------------------------------------------


/*
// Finds the smallest modulus for the local conditions with positive even valuation at the anisotropic primes.
mpz_class QuaternaryExceptions::_MakeNiceLocalModulus() const {
  
  // Declare some variables
  mpz_class new_mod = 1;
  vector< vector<long> > local_repn_vec;
  local_repn_vec = _local_conditions.Get_Local_Mod_Vector();
  vector<long> aniso_vec;
  aniso_vec = _local_conditions.Get_Anisotropic_Primes();
  
  
  // Find the smallest power of p needed in the local modulus
  for(long i=0; i<local_repn_vector.size(); i++) {
    long p_pow = 0;
    mpz_class p_part = 1;
    mpz_class p = local_repn_vec[i][0];
    
    if (p==2) {    // If p = 2
      p_pow = max(p_pow, local_repn_vec[i][1]);
      p_pow = max(p_pow, local_repn_vec[i][2]);
      p_pow = max(p_pow, local_repn_vec[i][3]);
      p_pow = max(p_pow, local_repn_vec[i][4]);

      p_pow = max(p_pow, local_repn_vec[i][5] + 1);
      p_pow = max(p_pow, local_repn_vec[i][6] + 1);
      p_pow = max(p_pow, local_repn_vec[i][7] + 1);
      p_pow = max(p_pow, local_repn_vec[i][8] + 1);

      p_pow = p_pow + 3;  // Since we start mod 8 = 2^3
    }
    else {        // If p > 2
      p_pow = max(p_pow, local_repn_vec[i][1]);
      p_pow = max(p_pow, local_repn_vec[i][2]);
      p_pow = max(p_pow, local_repn_vec[i][3] + 1);
      p_pow = max(p_pow, local_repn_vec[i][4] + 1);

      p_pow = p_pow + 1;  // Since we start mod p 

    }

    // Make the p-part of the modulus
    new_mod *= p ^ p_pow;  
  }


  // For anisotropic primes p, be sure p_pow is even and positive   
  for(long j=0; j < aniso_vec.size(); j++) {  
    unsigned long v = Valuation(new_mod, aniso_vec[j]);
    if (v == 0)
      new_mod = new_mod * aniso_vec[j] * aniso_vec[j];  // If p isn't there, introduce p^2.
    else if (v % 2 != 0)
      new_mod = new_mod * aniso_vec[j];                 // If p is there but an odd power, multiply by p.
  }



  // Return the nice modulus
  return new_mod;

}
*/  





// Finds the exceptions of a ternary form up to some bound

void QuaternaryExceptions::_QuaternaryExceptionsUpTo(const mpz_class & desired_precision, 
					       const string & ternary_exceptions_dir) {

  // For convenience, we define              ****** are these right???  don't we need half of a, b, and c??? ******
  mpz_class a = _T(1,1)/2;
  mpz_class b = _T(1,2);
  mpz_class c = _T(1,3);
  mpz_class d = _T(1,4);
  mpz_class e = _T(2,2)/2;
  mpz_class f = _T(2,3);
  mpz_class g = _T(2,4);
  mpz_class h = _T(3,3)/2;
  mpz_class i = _T(3,4);
  mpz_class j = _T(4,4)/2;



  mpz_class disc = 16*a*e*h*j - 4*a*e*i*i - 4*a*f*f*j + 4*a*f*g*i - 4*a*g*g*h - 4*b*b*h*j + b*b*i*i 
    + 4*b*f*c*j - 2*b*d*f*i - 2*b*g*c*i + 4*b*d*g*h - 4*e*c*c*j + 4*c*d*e*i + c*c*g*g - 2*c*d*g*f 
    - 4*d*d*e*h + d*d*f*f;

  mpz_class w11 = 8*e*h*j - 2*e*i*i - 2*f*f*j + 2*f*g*i - 2*g*g*h;
  mpz_class w22 = 8*a*h*j - 2*a*i*i - 2*c*d*j + 2*d*c*i - 2*d*d*h;
  mpz_class w33 = 8*a*e*j - 2*a*g*g - 2*b*b*j + 2*b*d*g - 2*d*d*e;
  mpz_class w44 = 8*a*e*h - 2*a*f*f - 2*b*b*h + 2*b*c*f - 2*c*c*e;


  /*
  // DIAGNOSTIC
  cout << " a = " << a << endl;
  cout << " b = " << b << endl;
  cout << " c = " << c << endl;
  cout << " d = " << d << endl;
  cout << " e = " << e << endl;
  cout << " f = " << f << endl;
  cout << " g = " << g << endl;
  cout << " h = " << h << endl;
  cout << " i = " << i << endl;
  cout << " j = " << j << endl;
  cout << endl;
  cout << " disc = " << disc << endl;
  cout << " w11 = " << w11 << endl;
  cout << " w22 = " << w22 << endl;
  cout << " w33 = " << w33 << endl;
  cout << " w44 = " << w44 << endl;
  cout << endl;
  */



  // Check the directory exists (If not, then forget it.)
  string Valid_Exceptions_dir;
  if (DirectoryExists(ternary_exceptions_dir) == true)
    Valid_Exceptions_dir = ternary_exceptions_dir;
  else {
    cout << " Warning: The directory passed to store ternary exceptions " << ternary_exceptions_dir << " is invalid."; 
    Valid_Exceptions_dir = "";
  }

  // Attempt to read previous results (if a valid directory is given)   
  // ----------------------------------------------------------------
  if (Valid_Exceptions_dir != "") {
 
   // Make the path to the Quaternary Exception file 
    string ExceptionFilename;
    ExceptionFilename = Valid_Exceptions_dir + "Quaternary_Exceptions__" + _T.QF_String() + ".txt";
    
  
    // Read the exceptions
    ifstream exception_file;
    exception_file.open(ExceptionFilename.c_str());
    if ((ternary_exceptions_dir != "") && (exception_file.is_open() == true)) {

      // No need to compute anything, so forget it.
      exception_file.close();
      
      // Read the series (which we know now exists!)
      cout << " Old _precision = " << _precision << endl;
      //      set<mpz_class> new_exception_set;
      //      mpz_class new_precision;
      _ReadExceptions(ExceptionFilename, _exception_set, _precision);      
      //      _exception_set = new_exception_set;
      //      _precision = new_precision;
      cout << " New _precision = " << _precision << endl;

    }
    
  }
  




  // Set variables for the range of numbers to check:
  // -------------------------------------------------  
  
  //  mpz_class Bunch_Size = 1000000;     // One million
  //  mpz_class Bunch_Size = 100000;      // One Hundred Thousand
  mpz_class Bunch_Size = 1000;         // 1 Thousand
  mpz_class Save_Increment = 100000;   // Save the computation after every hundred thousand numbers.
  mpz_class Report_Size_Drop = 10000;  // Report when the truant drops by 10 thousand

  mpz_class Checking_Min = _precision + 1;
  mpz_class Save_Target = Checking_Min + Save_Increment;  

  
  set<mpz_class> eligible_set;


  // Break our computation into bunches if necessary
  while (Checking_Min <= desired_precision) {

    mpz_class Checking_Max = min(desired_precision, Checking_Min + Bunch_Size);   // Says the check the range Checking_Min --> Checking_Max
    cout << endl;
    cout << " Checking from " << Checking_Min << " to " << Checking_Max << "  -----  " << Time();
    


    // Make a set of eligible numbers:
    // -------------------------------

    // /*  
    // DIAGNOSTIC
    cout << " Generating a batch of possible (locally represented) exceptions" << endl;
    // */
    
    // Approach #1:  (Lazy approach)
    for(mpz_class num = Checking_Min; num <= Checking_Max; num++) {
      if (_local_conditions.IsLocallyRepresented(num) == true)
	eligible_set.insert(num);
    }
    
    /*
    // Approach #2:  (Efficient approach)
    mpz_class nice_mod; 
    nice_mod = MakeNiceLocalModulus(local_conditions);  
    
    for(long mod_class = 0; mod_class < nice_mod; mod_class++) {      // This adds all locally represented numbers coming from this class
    
    // 
    
    
    }
    */
    
    
    // /*   
    // DIAGNOSTIC
    cout << " Finished generating a batch of possible (locally represented) exceptions from " << Checking_Min << " to " << Checking_Max << endl;
    cout << " There are " << eligible_set.size() << " of them." << endl;
    //cout << "   They are: " << endl << eligible_set << endl;
    // */

    
    
    
    // Compute the bounds on vectors to check:
    // ---------------------------------------
    cout << "   Checking possible (locally represented) exceptions from " << Checking_Min << " to " << Checking_Max << " ---- " << Time();

    
    // Approach #1: From Jagy
    mpz_class x_Max, y_Max, z_Max, w_Max;
    x_Max = ceil(sqrt(abs(Checking_Max.get_d() * w11.get_d() / disc.get_d())));
    y_Max = ceil(sqrt(abs(Checking_Max.get_d() * w22.get_d() / disc.get_d())));
    z_Max = ceil(sqrt(abs(Checking_Max.get_d() * w33.get_d() / disc.get_d())));
    w_Max = ceil(sqrt(abs(Checking_Max.get_d() * w44.get_d() / disc.get_d())));

    //    /*  
    // DIAGNOSTIC 
    cout << " w11 = " << w11 << "  " << w11.get_d() << endl;
    cout << " w22 = " << w22 << "  " << w22.get_d() << endl;
    cout << " w33 = " << w33 << "  " << w33.get_d() << endl;
    cout << " w44 = " << w44 << "  " << w44.get_d() << endl;
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
    cout << " w_Max = " << w_Max << endl;
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
    
    mpz_class WA, WB, WC;
    mpz_class Big_Disc, Small_Disc;
    mpz_class Sqrt_Big_Disc, Sqrt_Small_Disc;
    mpz_class w1, w2, w3, w4;
    
    mpz_class last_size = eligible_set.size();
    mpz_class diag, B, C, D, F, G, I;
    mpz_class I1, G1, D1, diag1;
    
    for(mpz_class x=0; (x <= x_Max) && (eligible_set.empty() == false); x++) 
      
      for(mpz_class y=0; (y <= y_Max) && (eligible_set.empty() == false); y++) 
	
	for(mpz_class z=z_Max; (z >= 0) && (eligible_set.empty() == false); z--) {
	  
	  // Computing the coefficients of the polynomial in w (given x,y,z)
	  WC = a*x*x + b*x*y + c*x*z + e*y*y + f*y*z + h*z*z;       // Temporary variable
	  WB = d*x + g*y + i*z;                                     // Temporary variable
	  WA = j;                                                   // Temporary variable
	  
	  
	  // Compute Big_Disc, Small_Disc, and TRUANT
	  Big_Disc = WB*WB  - 4*WA*(WC - MAX_TARGET);
	  Small_Disc = WB*WB  - 4*WA*(WC - TRUANT);
	  TRUANT = *eligible_set.begin();
	  
	  // Compute the bounds z1 -> z4
	  if (Big_Disc >= 0) {
	    Sqrt_Big_Disc = _Jagy_IntSqrt(Big_Disc);         // Temporary Variable
	    
	    w4 = ( -WB + Sqrt_Big_Disc ) / (2*WA);
	    w1 = ( -WB - Sqrt_Big_Disc ) / (2*WA);
	    
	    w3 = -WB / (2*WA);
	    w2 = w3;
	    
	    if (Small_Disc >= 0) {	
	      Sqrt_Small_Disc = _Jagy_IntSqrt(Small_Disc);   // Temporary Variable
	      w3 = ( -WB + Sqrt_Small_Disc ) / (2*WA);
	      w2 = ( -WB - Sqrt_Small_Disc ) / (2*WA);
	    }
	    
	  }


	  // Prepare for the 2 loops:
	  F = f*y*z;
	  B = b*x*y;
	  C = c*x*z;
	  I1 = i*z;
	  G1 = g*y;
	  D1 = d*x;
	  diag1 = a*x*x + e*y*y + h*z*z;

	  
	  // 2 Loops for w:  w3 --> w4  and  w1 --> w2 
	  // ------------------------------------------
	  for(mpz_class w=w3; (w <= w4) && (eligible_set.empty() == false); w++) {
	    
	    // Compute the diagonal and off-diagonal terms                           // <--- Can speed this up! =)
	    diag = diag1 + j*w*w;
	    I = I1*w;
	    G = G1*w;
	    D = D1*w;
	    
	    
	    
	    //	    /*	  
	    // DIAGNOSTIC:
	    mpz_class test_num(22);

	    if ((diag + B + C + D + F + G + I == test_num) 
	    || (diag - B - C - D + F + G + I == test_num) 
	    || (diag - B + C + D - F - G + I == test_num) 
	    || (diag + B - C + D - F + G - I == test_num) 
	    || (diag + B + C - D + F - G - I == test_num) 
	    || (diag + B - C - D - F - G + I == test_num) 
	    || (diag - B + C - D - F + G - I == test_num) 
	    || (diag - B - C + D + F - G - I == test_num)){
	    cout << endl;
	    cout << " Found the test number " << test_num << endl;
	    cout << "   diag + B + C + D + F + G + I = " << (diag + B + C + D + F + G + I) << endl;
	    cout << "   diag - B - C - D + F + G + I = " << (diag - B - C - D + F + G + I) << endl;
	    cout << "   diag - B + C + D - F - G + I = " << (diag - B + C + D - F - G + I) << endl;
	    cout << "   diag + B - C + D - F + G - I = " << (diag + B - C + D - F + G - I) << endl;
	    cout << "   diag + B + C - D + F - G - I = " << (diag + B + C - D + F - G - I) << endl;
	    cout << "   diag + B - C - D - F - G + I = " << (diag + B - C - D - F - G + I) << endl;
	    cout << "   diag - B + C - D - F + G - I = " << (diag - B + C - D - F + G - I) << endl;
	    cout << "   diag - B - C + D + F - G - I = " << (diag - B - C + D + F - G - I) << endl;
	    cout << endl;	      
	    cout << "   x = " << x << "   y = " << y << "   z = " << z << "   w = " << w << endl;
	    cout << "   a = " << a << "   b = " << b << "   c = " << c << "   d = " << d << "   e = " << e 
	    << "   f = " << f << "   g = " << g << "   h = " << h << "   i = " << i << "   j = " << j << endl; 
	    cout << "   diag = " << diag << "   B = " << B << "   C = " << C << "   D = " << D 
	    << "   F = " << F << "   G = " << G << "   I = " << I << endl;
	    cout << endl;	      	      
	    }
	    //	    */
	    
	    
	    
	    // Compute and remove the 8 values
	    eligible_set.erase(diag + B + C + D + F + G + I);  // No switching
	    
	    eligible_set.erase(diag - B - C - D + F + G + I);  // switch x -> -x
	    eligible_set.erase(diag - B + C + D - F - G + I);  // switch y -> -y
	    eligible_set.erase(diag + B - C + D - F + G - I);  // switch z -> -z
	    eligible_set.erase(diag + B + C - D + F - G - I);  // switch w -> -w
	    
	    eligible_set.erase(diag + B - C - D - F - G + I);  // switch (x,y) -> -(x,y)
	    eligible_set.erase(diag - B + C - D - F + G - I);  // switch (x,z) -> -(x,z)
	    eligible_set.erase(diag - B - C + D + F - G - I);  // switch (x,w) -> -(x,w)
	    
	    
	    // Check if the size has dropped by Report_Size_Drop (default is 10,000)
	    if (eligible_set.size() < last_size - Report_Size_Drop) {
	      cout << "  " << Time() <<" --  (1) The current eligible exception set has size " << eligible_set.size() << endl;
	      last_size = eligible_set.size();
	    }
	  }
	  
	  
	  for(mpz_class w=w1; (w <= w2) && (eligible_set.empty() == false); w++) {
	    
	    
	    // Compute the diagonal and off-diagonal terms                           // <--- Can speed this up! =)
	    diag = diag1 + j*w*w;
	    I = I1*w;
	    G = G1*w;
	    D = D1*w;
	    
	    
	    // Compute and remove the 8 values
	    eligible_set.erase(diag + B + C + D + F + G + I);  // No switching
	    
	    eligible_set.erase(diag - B - C - D + F + G + I);  // switch x -> -x
	    eligible_set.erase(diag - B + C + D - F - G + I);  // switch y -> -y
	    eligible_set.erase(diag + B - C + D - F + G - I);  // switch z -> -z
	    eligible_set.erase(diag + B + C - D + F - G - I);  // switch w -> -w
	    
	    eligible_set.erase(diag + B - C - D - F - G + I);  // switch (x,y) -> -(x,y)
	    eligible_set.erase(diag - B + C - D - F + G - I);  // switch (x,z) -> -(x,z)
	    eligible_set.erase(diag - B - C + D + F - G - I);  // switch (x,w) -> -(x,w)
	    
	    
	    // Check if the size has dropped by Report_Size_Drop (default is 10,000)
	    if (eligible_set.size() < last_size - Report_Size_Drop) {
	      cout << "  " << Time() <<" --  (2) The current eligible exception set has size " << eligible_set.size() << endl;
	      last_size = eligible_set.size();
	    }
	    
	    
	  }	
	}
    
    
    
    
    // Store the cumulative results in _exception_set and _precision
    _exception_set.insert(eligible_set.begin(), eligible_set.end());
    _precision = Checking_Max;
    
    
    // Write the cumulative resuts (if a valid directory was given) if we're big enough or done.
    // -----------------------------------------------------------------------------------------
    
    if ((Checking_Max >= desired_precision) || (Checking_Max >= Save_Target)) {

      // Make the path to the Quaternary Exception file 
      string ExceptionFilename;
      ExceptionFilename = Valid_Exceptions_dir + "Quaternary_Exceptions__" + _T.QF_String() + ".txt";
      
      // Write the results
      _WriteExceptions(ExceptionFilename, _exception_set, _precision);      

      // Increment the target for the next save operation
      Save_Target += Save_Increment;

    }


    // Increment the bottom of the checking range
    Checking_Min += Bunch_Size;

    
  }   // End of checking in bunches -- All numbers checked! =)
  
  
  
}





// Checks if a given number is represented by a quaternary form
bool QuaternaryExceptions::_IsQuaternaryException(const mpz_class & num) {

  // For convenience, we define              ****** are these right???  don't we need half of a, b, and c??? ******
  mpz_class a = _T(1,1)/2;
  mpz_class b = _T(1,2);
  mpz_class c = _T(1,3);
  mpz_class d = _T(1,4);
  mpz_class e = _T(2,2)/2;
  mpz_class f = _T(2,3);
  mpz_class g = _T(2,4);
  mpz_class h = _T(3,3)/2;
  mpz_class i = _T(3,4);
  mpz_class j = _T(4,4)/2;



  mpz_class disc = 16*a*e*h*j - 4*a*e*i*i - 4*a*f*f*j + 4*a*f*g*i - 4*a*g*g*h - 4*b*b*h*j + b*b*i*i 
    + 4*b*f*c*j - 2*b*d*f*i - 2*b*g*c*i + 4*b*d*g*h - 4*e*c*c*j + 4*c*d*e*i + c*c*g*g - 2*c*d*g*f 
    - 4*d*d*e*h + d*d*f*f;

  mpz_class w11 = 8*e*h*j - 2*e*i*i - 2*f*f*j + 2*f*g*i - 2*g*g*h;
  mpz_class w22 = 8*a*h*j - 2*a*i*i - 2*c*c*j + 2*d*c*i - 2*d*d*h;
  mpz_class w33 = 8*a*e*j - 2*a*g*g - 2*b*b*j + 2*b*d*g - 2*d*d*e;
  mpz_class w44 = 8*a*e*h - 2*a*f*f - 2*b*b*h + 2*b*c*f - 2*c*c*e;


  /*
  // DIAGNOSTIC
  cout << " a = " << a << endl;
  cout << " b = " << b << endl;
  cout << " c = " << c << endl;
  cout << " d = " << d << endl;
  cout << " e = " << e << endl;
  cout << " f = " << f << endl;
  cout << " g = " << g << endl;
  cout << " h = " << h << endl;
  cout << " i = " << i << endl;
  cout << " j = " << j << endl;
  cout << endl;
  cout << " disc = " << disc << endl;
  cout << " w11 = " << w11 << endl;
  cout << " w22 = " << w22 << endl;
  cout << " w33 = " << w33 << endl;
  cout << " w44 = " << w44 << endl;
  cout << endl;
  */





  // Set variables for the range of numbers to check:
  // -------------------------------------------------    
  set<mpz_class> eligible_set;
  
  
  
  // Make a set of eligible numbers:
  // -------------------------------
  
  
  // Put the number num in the eligible set
  if (_local_conditions.IsLocallyRepresented(num) == true)
    eligible_set.insert(num);
  
  
  
  // /*   
  // DIAGNOSTIC
  cout << "   Checking the set of numbers: " << eligible_set << endl;
  // */
  
  
  
  
  // Compute the bounds on vectors to check:
  // ---------------------------------------
  
  
  // Approach #1: From Jagy
  mpz_class x_Max, y_Max, z_Max, w_Max;
  x_Max = ceil(sqrt(abs(2*num.get_d() * w11.get_d() / disc.get_d())));
  y_Max = ceil(sqrt(abs(2*num.get_d() * w22.get_d() / disc.get_d())));
  z_Max = ceil(sqrt(abs(2*num.get_d() * w33.get_d() / disc.get_d())));
  w_Max = ceil(sqrt(abs(2*num.get_d() * w44.get_d() / disc.get_d())));
  
  // /*  
  // DIAGNOSTIC 
  if (num == 406) {
    cout << " w11 = " << w11 << "  " << w11.get_d() << endl;
    cout << " w22 = " << w22 << "  " << w22.get_d() << endl;
    cout << " w33 = " << w33 << "  " << w33.get_d() << endl;
    cout << " w44 = " << w44 << "  " << w44.get_d() << endl;
    cout << " disc = " << disc << "  " << disc.get_d() << endl;
    cout << endl;
    cout << " Using the bounds: " << endl;
    cout << " x_Max = " << x_Max << endl;
    cout << " y_Max = " << y_Max << endl;
    cout << " z_Max = " << z_Max << endl;
    cout << " w_Max = " << w_Max << endl;
    cout << endl;  
  }
  //  */
  
  
  
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
  mpz_class TRUANT, MAX_TARGET = num;                       // <--- NEED TO SET THESE! =)
  
  mpz_class WA, WB, WC, W_Disc;
  mpz_class Big_Disc, Small_Disc;
  mpz_class Sqrt_Big_Disc, Sqrt_Small_Disc;
  mpz_class w1, w2, w3, w4;
  
  mpz_class last_size = eligible_set.size();
  mpz_class diag, B, C, D, F, G, I;
  mpz_class I1, G1, D1, diag1;
  
  for(mpz_class x=0; (x <= x_Max) && (eligible_set.empty() == false); x++) 
    
    for(mpz_class y=0; (y <= y_Max) && (eligible_set.empty() == false); y++) 
      
      for(mpz_class z=z_Max; (z >= 0) && (eligible_set.empty() == false); z--) {

	// Loop through all possible signs of (x,y,z) since z is always positive! =)
	mpz_class X,Y,Z;
	X = x;  Y = y;  Z = z;
	mpz_class w;

	//	for(long sign_x=0; sign_x <=1; sign_x++) {
	//	  X *= -1;
	for(long sign_y=0; sign_y <=1; sign_y++) {
	  Y *= -1;
	  for(long sign_z=0; sign_z <=1; sign_z++) {
	    Z *= -1;
	    
	    for(long sign_w=0; sign_w <=1; sign_w++) {
	      
	      
	      
	      // Computing the coefficients of the polynomial in w (given x,y,z)
	      WC = a*X*X + b*X*Y + c*X*Z + e*Y*Y + f*Y*Z + h*Z*Z;       // Temporary variable
	      WB = d*X + g*Y + i*Z;                                     // Temporary variable
	      WA = j;                                                   // Temporary variable
	      
	      // Solve for w, and check the solution.
	      W_Disc = WB*WB  - 4*WA*(WC - num);
	      if (W_Disc >= 0) { 
		if (sign_w == 0)
		  w = ( -WB + sqrt(W_Disc) ) / (2*WA);
		else
		  w = ( -WB - sqrt(W_Disc) ) / (2*WA);
	      }

	      /*
	      // DIAGNOSTIC for the missing representation (17,6,5,0) of 406 for form # 6230:
	      if ((num == 406) && (x == 17) && (y == 6) && (z==5)) {
		cout << endl << endl;
		cout << " For the missing representation (17,6,5,0), we are using the bounds: " << endl;
		cout << " w1 = " << w1 << endl;
		cout << " w2 = " << w2 << endl;
		cout << " w3 = " << w3 << endl;
		cout << " w4 = " << w4 << endl;
		cout << endl;
		cout << " WA = " << WA << endl;
		cout << " WB = " << WB << endl;
		cout << " WC = " << WC << endl;
		cout << endl;
		cout << " W_Disc = " << W_Disc << endl;
		cout << " sqrt(W_Disc) = " << sqrt(W_Disc) << endl;
		cout << " -WB + sqrt(W_Disc) = " << (-WB + sqrt(W_Disc)) << endl;
		cout << " -WB - sqrt(W_Disc) = " << (-WB - sqrt(W_Disc)) << endl;
		cout << " (-WB + sqrt(W_Disc)) / (2*WA) = " << ((-WB + sqrt(W_Disc)) / (2*WA)) << endl;
		cout << " (-WB - sqrt(W_Disc)) / (2*WA) = " << ((-WB - sqrt(W_Disc)) / (2*WA)) << endl;
		cout << endl;
		cout << " X = " << X << endl;
		cout << " Y = " << Y << endl;
		cout << " Z = " << Z << endl;
		cout << " w = " << w << endl;
		cout << endl;
	      }
	      */


	      
	      // Remove the associated value
	      eligible_set.erase(WA*w*w + WB*w + WC);       // IT MIGHT BE FASTER TO RETURN DIRECTLY IF THE VALUE num IS ATTAINED! =)
	      
	    }
	  }
	}
	
      }


	
  
  /*
  // DIAGNOSTIC:
  if (num == 406) 
    cout << "EXITING WITH eligible_set.empty() = " << eligible_set.empty() << endl;
  */  

  
  // Return the result
  if (eligible_set.empty() == true)
    return false;
  else
    return true;
  
}





// Jagy's IntSqrt routine used for his ternary exception finder. =)
// (This returns the largest integer x so that x^2 < n.)                    <=====  CHECK THIS... the end is confusing...
inline 
mpz_class QuaternaryExceptions::_Jagy_IntSqrt(const mpz_class & n) const {
  
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








//////////////////////////////////
// Reads ternary exceptions from a file
//////////////////////////////////
void QuaternaryExceptions::_ReadExceptions(const string & ExceptionFilename, set<mpz_class> & new_exception_set, mpz_class & new_precision) const {



  cout << " Reading the file: " << ExceptionFilename <<  endl;


  // Try to open the file
  ifstream exceptionfile;
  exceptionfile.open(ExceptionFilename.c_str(), ios::in);

  // Abort if we fail... =(
  if (! exceptionfile.is_open())
    { cout << "ReadExceptions Error: Error opening file"; exit(1); }


  // Read the exceptions
  char c;
  mpz_class file_precision;
  set<mpz_class> file_exception_set;

  exceptionfile >> file_precision;               // Read the precision
  //  file_exception_set = ReadSet_mpz_class(exceptionfile);   // Read the set of exceptions
  exceptionfile >> file_exception_set;           // Read the set of exceptions

  // Sanity Checks
  assert( file_precision >= 0 );                      // Check the precision isn't negative
                                                                              // How to Check the biggest exception is less than the precision ???
  cout << " file_precision = " << file_precision << endl;
  cout << " file_exception_set = " << file_exception_set << endl;


  // Close the file
  exceptionfile.close();


  // Return the precision and exception set
  new_precision = file_precision;
  new_exception_set = file_exception_set;

}                








////////////////////////////////
// Writes the ternary exceptions to a file
////////////////////////////////
void QuaternaryExceptions::_WriteExceptions(const string & ExceptionFilename, const set<mpz_class> & new_exception_set, const mpz_class & new_precision) const {

  // Try to open the file
  ofstream exceptionfile;
  exceptionfile.open(ExceptionFilename.c_str(), ios::out);

  // Abort if we fail... =(
  if (! exceptionfile.is_open())
    { cout << "WriteExceptions Error: Error opening file"; exit (1); }

  // Write the file 
  exceptionfile << _precision << endl;
  exceptionfile << _exception_set << endl;

  // Close the file
  exceptionfile.close();

}

