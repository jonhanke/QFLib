


////////////////////////////////////////////////////////////////
// Writes the boolean ternary theta function as a binary file //
// NOTE: This is just a copy of the routine in maketheta.cc!  //
// (with some minor modifications to the directory and filename...)
////////////////////////////////////////////////////////////////

void boolean_ternary_theta::_WriteTernaryThetaBinary_client(char* hostname, long slice) {

  // *** WARNING:  We assume that theta_bin has size precision/32 
  // ***           (which could also be coded as precision >> 5).

  // Notation:
  // ---------
  // precision of the theta function ==> ... + O(x^(precision + 1))
  //
  // Q = [ a, b, c ]
  //     [ 0, d, e ]
  //     [ 0, 0, f ]

  //  extern char DIST_DIR[];  // This is the global directory for the temporary distributed theta function files
  
  // Make the filename "Slice_#__HOST__ternary_bool_theta__a_b_c_d_e_f__PRECISION.bindata"
  char filename[200];  
  sprintf(filename, "%sSlice_%d__%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, slice, hostname, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 

  // Make the filename "Slice_#__HOST__ternary_bool_theta__a_b_c_d_e_f__PRECISION.bindata"
  char temp_filename[200];  
  sprintf(temp_filename, "%sTEMP__Slice_%d__%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, slice, hostname, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 


  // Write the temporary file
  write(temp_filename);

  // Move the temporary file to it's real location
  char mv_command[600];  
  sprintf(mv_command, "mv %s %s", temp_filename, filename); 
  system(mv_command);

}





// Main program (server) for computing the distributed theta function

void boolean_ternary_theta::compute_client(const long slice_number) {

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

  
  // Find the Cholesky decomposions
  double Cholesky_temp[6];
  CholeskyTernary(QQ, Cholesky_temp);



  // Get the hostname
  char hostname[255];
  int returncode;
  returncode = gethostname(hostname, 255);  


  // Start the computation for this slice
  cout << "\n\nStarting the slice" << endl << endl;
  _compute_client_slice(slice_number);
  cout << "\n\nFinishing the slice" << endl << endl;


  // DIAGNOSTIC -- to check the computed values before writing
  for(long kk=0; kk < precision(); kk++)
    cout << " kk =  " << kk << "    value = " << get_value(kk) << endl;


  // Write the output to the Distributed directory
  _WriteTernaryThetaBinary_client(hostname, slice_number); 

}













// ==============================================================================================









////////////////////////////////////////////////////////////////////////////////////////
// Finds all vectors with Q(x) <= C, so long as the slice_size is less than slice_min //
////////////////////////////////////////////////////////////////////////////////////////

long boolean_ternary_theta::_compute_smallslices(const unsigned long SLICE_MIN) { 
  
  double Cholesky_temp[6];
  CholeskyTernary(QQ, Cholesky_temp);

  vector<long> empty_vec;

  return _compute_with_slices(Cholesky_temp, SLICE_MIN, empty_vec);
  
}




//////////////////////////////////////////////////////
// Find all vectors with Q(x) <= C in a given slice //
//////////////////////////////////////////////////////

void boolean_ternary_theta::_compute_client_slice(const long slice_num) { 
  
  double Cholesky_temp[6];
  CholeskyTernary(QQ, Cholesky_temp);
  vector<long> slice_vec(1, slice_num);

  _compute_with_slices(Cholesky_temp, 0, slice_vec);
  
}





////////////////////////////////////////////////////////////////////////////////////////
// Finds all vectors with Q(x) <= C, 
// so long as the slice_size is less than slice_min //
// or the 
////////////////////////////////////////////////////////////////////////////////////////

long boolean_ternary_theta::_compute_with_slices(const double Cholesky[], const unsigned long SLICE_MIN, const vector<long> slice_vec) {
  //long FastBinaryThetaInitial(unsigned long theta[], double QQ[], double C, long slice_min) { 

  if ((SLICE_MIN != 0) && (slice_vec.size() != 0)) {
    cout << "ERROR IN _compute_with_slices:" 
	 << " Can't decide if we're doing one slice, or all small slices..." << endl;  
    abort();
  }

  // Check all slices are <= zero
  for (long k=0; k < slice_vec.size(); k++)
    if (slice_vec[k] > 1) {
      cout << "ERROR IN _compute_with_slices:" 
	   << " The slices to compute must be <= 0 (or = 1 for the _compute_smallslices front-end)...." << endl;  
      abort();
    }


  cout << endl << "Entering _compute_with_slices(...)" << endl;
  //  cout << " Using C = " << C << endl;


  // ERROR CHECKING: Check that C+1 fits in an unsigned long.



  // Zero out the theta function
  // for (unsigned long i=0; i < precision(); i++)
  // --WGM
  for (unsigned long i=0; i < _length(); i++)
    _theta[i] = 0;


  const long n = 3;

  // Make the (constant) matrix Q from Cholesky  -- NOTE: We only use indices 1--3.
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
  T[i] = (double) precision();
  U[i] = 0;
  double Z;
  vector<long> L(n+1, 0);
  vector<long> x(n+1, 0);
  unsigned long slice_counter = 0;   // Counts the number of vectors in a given slice

  cout << "Finished the intialization." << endl;


  // 2. Compute bounds
  Z = sqrt(T[i] / Q[i][i]);
  L[i] = long(floor(Z - U[i]));  // Check this is ok...
  /*
  cout << " L[i] float gives :     " << floor(Z - U[i].get_d()) << endl;
  cout << " L[i] mpz_class gives : " << L[i] << endl;
  */
  
  if (slice_vec.size() != 0)
    x[i] = slice_vec[0] - 1;    // This initializes x[n] to be the desired slice, we we're doing one slice
                                // This allows for the increment in the beginning of the loop (in 3a).
  else
    x[i] = long(ceil(-Z - U[i]) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i].get_d() -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */

  cout << "Finished computing the initial bounds." << endl;
  
  bool done_flag = false;
  bool slice_overflow_flag = false;
  long slice_ptr = 0;                 // Tells which is the current slice being counted (in the client routine)

  double Q_val_double;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...
  cout << "Finished declaring some variables." << endl;
  
  
  // Big loop which runs through all vectors
  // (loop if we're not done, and either we're on the last slice or slice_counter is big enough to distribute the computation)
  // (so this doesn't finish the slice 
  while ((done_flag == false) && (slice_overflow_flag == false)) {
    
    
    // Loop through until we get to i=1 (so we define a vector x)
    do {
      
      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }


      // When i==n, we're in a different slice!
      if (i == n) {
	
	// If we're counting the size of the slices, then clear the slice_counter
	if (SLICE_MIN > 0) {
	  cout << "slice_counter = " << slice_counter << endl;
	  if (slice_counter >= SLICE_MIN)
	    slice_overflow_flag = true;
	  
	  slice_counter = 0;
	}
	
	// ----------------------------------
	
	// If we're only checking a few slices and we overflow the last slice,
	// then move to the next slice (so we force 3a) unless we're done.
	if ((slice_vec.size() != 0) && (x[n] != slice_vec[slice_ptr])) {
	  if (slice_ptr != (slice_vec.size() - 1)) {
	    // /*
	    cout << " aaa " << endl;
	    cout << " slice_ptr = " << slice_ptr << endl;
	    cout << " slice_vec.size() = " << slice_vec.size() << endl;
	    cout << " x[n] = " << x[n] << "  slice_vec[slice_ptr] = " << slice_vec[slice_ptr] << endl;
	    // */
	    slice_ptr++;
	    x[n] = slice_vec[slice_ptr];
	  }
	  else {
	    //	    /*
	    cout << " bbb " << endl;
	    cout << " slice_ptr = " << slice_ptr << endl;
	    cout << " slice_vec.size() = " << slice_vec.size() << endl;
	    cout << " x[n] = " << x[n] << "  slice_vec[slice_ptr] = " << slice_vec[slice_ptr] << endl;
	    //	    */
	    done_flag = true;
	  }
	  
	}      
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
      
    } while ((i > 1) && (done_flag == false));
    
    
    // 4. Solution found (This happens when i=1)
    /*
    cout << " i = " << i << endl;
    cout << " done_flag = " << done_flag << endl;
    cout << " precision() = " << precision() << endl;
    cout << " x = [ " << x[1] << ", " << x[2] << ", " << x[3] << " ]" << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */
    Q_val_double = precision() - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    Q_val = (unsigned long) round(Q_val_double);

    //    cout << " Float = " << Q_val_double << "   Long = " << Q_val << "  XX " << endl;
    /*
    cout << " The float value is " << Q_val_double << endl;
    cout << " The associated long value is " << Q_val << endl;
    cout << endl;
    */



    /*
    cout << " i = " << i << endl;
    cout << " done_flag = " << done_flag << endl;
    cout << " precision() = " << precision() << endl;
    cout << " (Q_val <= precision()) = " << (Q_val <= precision()) << endl;
    */
    /*
    cout << " x = [ " << x[1] << ", " << x[2] << ", " << x[3] << " ]" << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */

    // Note: The last condition checks we haven't overflowed the current slice!
    if ((Q_val <= precision()) && (done_flag == false)) {
      //  cout << "Set the bit " << Q_val << endl << endl;
      set_value(Q_val);                                             // Set the appropriate bit
        //      theta[Q_val] = theta[Q_val] + 2;
      slice_counter++;                                             // Increment the slice_counter
    }
    

 
    // 5. Check if x = 0, for exit condition. =)
    long j=1;
    done_flag = true;
    while (j<=n) {
      if (x[j] != 0)
	done_flag = false;
      j++;
    }    
    
    // 5a. If we overflow the last allowed slice, then leave also.
    if ((slice_vec.size() != 0) && (x[n] != slice_vec[slice_ptr]) && (slice_ptr == (slice_vec.size() - 1)))
      done_flag = true;



  }
	 
	 
  cout << "Leaving FastBinaryThetaInitial" << endl << endl;


  if (slice_overflow_flag == true)
    return x[n];   // Return the value of the slice we need to compute next!
  else
    return 1;      // same as above, but if we're done, then the value of the next slice should be 1.
                   // (even though we don't need to compute it!)
}
