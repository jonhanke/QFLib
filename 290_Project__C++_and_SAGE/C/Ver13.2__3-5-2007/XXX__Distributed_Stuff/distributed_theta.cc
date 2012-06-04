

////////////////////////////////////////////////////////////////////////////////////////
// Finds all vectors with Q(x) <= C, so long as the slice_size is less than slice_min //
////////////////////////////////////////////////////////////////////////////////////////

long FastBinaryThetaInitial(unsigned long theta[], double QQ[], double C, long slice_min) { 



  // *** WARNING:  We assume that theta has size ceil(C / 32), so it can hold the precision. 

  

  cout << endl << "Entering FastBinaryThetaInitial" << endl;
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
  long slice_counter = 0;   // Counts the number of vectors in a given slice


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
  bool slice_overflow_flag = false;
  double Q_val_double;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...
  
  
  // Big loop which runs through all vectors
  // (loop if we're not done, and either we're on the last slice or slice_counter is big enough to distribute the computation)
  // (so this doesn't finish the slice 
  while ((done_flag == false) && (slice_overflow_flag == false)) {
    
    
    // Loop through until we get to i=1 (so we defined a vector x)
    do {
      
      // Check to see if we should clear the slice_counter 
      // (i.e., when i=1 since the increment is in the next section)
      if (i == n) {
	if (slice_counter >= slice_min)
	  slice_overflow_flag = true;

	slice_counter = 0;
      }      

      
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
      slice_counter ++;                                             // Increment the slice_counter
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


  if (slice_overflow_flag == true)
    return x[n];   // Return the value of the slice we need to compute next!
  else
    return 1;      // same as above, but if we're done, then the value of the next slice should be 1.
                   // (even though we don't need to compute it!)
}







///////////////////////////////////////////////////////////////////////////////////////
// Reads the boolean ternary theta function for a slice from the temporary directory //
///////////////////////////////////////////////////////////////////////////////////////

bool ReadTernaryThetaBinary_distributed_host(unsigned long theta_read[], long form[6], unsigned long precision, const char* host) {

  // *** WARNING:  We assume that theta_bin has size precision/32 
  // ***           (which could also be coded as precision >> 5).

  // Notation:
  // ---------
  // precision of the theta function ==> ... + O(x^(precision + 1))
  //
  // Q = [ a, b, c ]
  //     [ 0, d, e ]
  //     [ 0, 0, f ]


  extern char DIST_DIR[];  // This is the global directory for the temporary distributed theta function files
  
  
  // Make the filename for each host "HOST__ternary_bool_theta__a_b_c_d_e_f__precision.bindata"
  char filename[200];  
  sprintf(filename, "%s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, host, form[0], form[1], form[2], form[3], form[4], form[5], precision); 


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
  bool file_opened = filein.is_open();

  if (file_opened == true) {

    // Read the array
    filein.read( (char*) theta_read, sizeof(long) * n);
    
    // Close the file
    filein.close();
  } 

  return file_opened;

}





















// Main program (server) for computing the distributed theta function

void MakeTernaryThetaDistributed__by_x_n(long QQ[6], unsigned long Ternary_Precision) {

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
  
  // Distribute the computation after the slice size reaches this (now it's 10^7)
    long SLICE_MIN = 10000000;
  //  long SLICE_MIN = 1000;


    /*
  // Use the Boolean_ternary_theta methods to do this
  boolean_ternary_theta theta3(QQ, Ternary_Precision);

  cout << " Computing the theta function " << endl;
  PrintTime();

  theta3.compute();
  theta3.write();
    */


  double Cholesky_temp[6];
  CholeskyTernary(QQ, Cholesky_temp);
  

  cout << "Allocating " << ((Ternary_Precision >> 5) +1) << " unsigned longs." << endl;
  //  unsigned long theta[(Ternary_Precision >> 5) +1];  // <-- This was the old allocation procedure... fails for 10^7.
  //  unsigned long *theta = (unsigned long *) calloc(sizeof(long), 10000000);
  unsigned long *theta = (unsigned long *) calloc(sizeof(long), ((Ternary_Precision >> 5) +1));
  cout << " Finished allocating the space!" << endl;

  // Start the computation for small slices
  long slice_todo;
  slice_todo = FastBinaryThetaInitial(theta, Cholesky_temp, Ternary_Precision, SLICE_MIN);


  // DIAGNOSTIC:
  cout << " The last two longs are: " << endl;
  cout << "  i = " << ((Ternary_Precision >> 5) - 1) << " theta[i] = " << theta[(Ternary_Precision >> 5) - 1] << endl;
  cout << "  i = " << ((Ternary_Precision >> 5) + 0) << " theta[i] = " << theta[(Ternary_Precision >> 5) + 0] << endl;
  cout << endl;


  // Check if we need to distribute
  if (slice_todo < 1) {

    // Allocate a temporary theta vector for the incoming slice information
    cout << "Allocating another " << ((Ternary_Precision >> 5) +1) << " unsigned longs." << endl;
    unsigned long *theta_temp = (unsigned long *) calloc(sizeof(long), ((Ternary_Precision >> 5) +1));
    cout << " Finished allocating the space!" << endl;

        
    // Make a list of the available computers
    string hosts[] = {"tux1", "tux2", "tux3"};

    cout << " There are " << sizeof(hosts)/sizeof(string) << " hosts." << endl;
    long num_of_hosts = sizeof(hosts)/sizeof(string);

    //    char hosts[num_of_hosts][100] = 

        // Some other possibilities:
        // -------------------------
        //     vector<char[100]> hosts(num_of_hosts)= {"tux1", "tux2", "tux3"};
        //
        // vector<string> hosts(num_of_hosts);
        // hosts = {"tux1", "tux2", "tux3"};



    // Make a directory to look for the partial results
    extern char DIST_DIR[];  // This is the global directory for the temporary distributed theta function files


    // Determine the way of slicing into chunks
    long current_host_slice[num_of_hosts];  // Stores the slice currently being worked on for each host.
    char* host_assignment_time[num_of_hosts];  // Stores the time each host was assigned its slice.


    // Delegate chunks in order
    for(long i=0; i < num_of_hosts; i++) {

      // Make the delegation command "ssh HOST 'nice ternary_theta_client a b c d e f precision slice_num'"
      char spawn_command[300];
      sprintf(spawn_command, "ssh %s 'nice ternary_theta_client %d %d %d %d %d %d %d %d'",
	  hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision, slice_todo); 
      
      // Spawn the slice computation
      int return_code;
      return_code = system(spawn_command);


      // Record the information for that host
      current_host_slice[i] = slice_todo;
      time_t rawtime;
      time(&rawtime);
      host_assignment_time[i] = ctime(&rawtime);
      slice_todo++;

      cout << " Sent slice " << current_host_slice[i] << " to " << hosts[i] << " at " << host_assignment_time[i] << endl; 
      
    }

     

    // Look for results and re-delegate
    for(long i=0; i < num_of_hosts; i++) {

      if (ReadTernaryThetaBinary_distributed_host(theta_temp, QQ, Ternary_Precision, hosts[i].c_str()) == true) {

	cout << " Recovered the data from host " << hosts[i] << endl;

	// Combine the new series with the old one
	for (unsigned long j = 0; j < ((Ternary_Precision >> 5) +1); j++)
	  theta[j] = theta[j] | theta_temp[j];

	// Remove the temporary computation for that host
	char rm_filename[300];  
	sprintf(rm_filename, "rm -f %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
		DIST_DIR, hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision); 
	system("rm_filename");
	
	// Send another slice to that host             TODO: Check that there are still slices to send!!! =)
	char spawn_command[300];
	sprintf(spawn_command, "ssh %s 'nice ternary_theta_client %d %d %d %d %d %d %d %d'",
		hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision, slice_todo); 
	int return_code;
	return_code = system(spawn_command);

	// Record the information for that host
	current_host_slice[i] = slice_todo;
	time_t rawtime;
	time(&rawtime);
	host_assignment_time[i] = ctime(&rawtime);
	slice_todo++;

	cout << " Sent slice " << current_host_slice[i] << " to " << hosts[i] << " at " << host_assignment_time[i] << endl; 


	////////////////////////////////////////////////////////////////////////
	// Save the current state in two temporary files (theta and state file)
	////////////////////////////////////////////////////////////////////////
	WriteTernaryThetaBinary_client(theta, QQ, Ternary_Precision, "Server_TEMPORARY"); 

	// Make the filename "Server_TEMPORARY__ternary_bool_theta__a_b_c_d_e_f__precision.sliceinfo"
	char filename[200];  
	sprintf(filename, "%s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision); 
	

	// Open the file for writing and check it opened correctly
	ofstream fileout;
	fileout.open(filename, ios::out | ios::trunc);
	
	if (! fileout.is_open())
	  { cout << "Error opening output file " << filename; exit (1); }

	// Write the table of slices and times.
	for (long j=0; j < num_of_hosts ; j++) 
	  fileout << hosts[j] << "  /  " 
		  << current_host_slice[j] << "  /  " 
		  << host_assignment_time[j] << "  /  " << endl;
	

	// Close the file
	fileout.close();



	// Replace the older files with the newer ones (use mv)
	char mv_theta[600];  
	sprintf(mv_theta, "mv -f %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision,
		DIST_DIR, "Server_CURRENT", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision); 
	system(mv_theta);
	char mv_sliceinfo[600];  
	sprintf(mv_sliceinfo, "mv -f %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision,
		DIST_DIR, "Server_CURRENT", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], Ternary_Precision); 
	system(mv_sliceinfo);

      }

    }


    
    // Look for timeouts and re-delegate -- Still TO DO
    
    
    
  }
  else if (slice_todo > 1) 
    cout << " ERROR in MakeTernaryThetaDistributed__by_x_n: The current slice_todo " << slice_todo << " is > 1!" << endl;
  

  // Write the final output when we're done! =)
  WriteTernaryThetaBinary(theta, QQ, Ternary_Precision); 
}




