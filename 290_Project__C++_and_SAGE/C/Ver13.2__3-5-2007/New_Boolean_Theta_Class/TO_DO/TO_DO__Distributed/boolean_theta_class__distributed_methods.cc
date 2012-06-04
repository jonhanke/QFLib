

/////////////////////////////
// Checks if a file exists //
/////////////////////////////

bool FileExists(const char* filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}



// ------------------------------------------------------------



////////////////////////////////////////////////////////////////////////////////////////
// Finds all vectors with Q(x) <= C, so long as the slice_size is less than slice_min //
////////////////////////////////////////////////////////////////////////////////////////

long boolean_ternary_theta::_compute_smallslices(double Cholesky[], unsigned long SLICE_MIN) {
  //long FastBinaryThetaInitial(unsigned long theta[], double QQ[], double C, long slice_min) { 


  cout << endl << "Entering FastBinaryThetaInitial" << endl;
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

  x[i] = long(ceil(-Z - U[i]) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i].get_d() -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */

  cout << "Finished computing the initial bounds." << endl;
  
  bool done_flag = false;
  bool slice_overflow_flag = false;
  double Q_val_double;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...
  cout << "Finished declaring some variables." << endl;
  
  
  // Big loop which runs through all vectors
  // (loop if we're not done, and either we're on the last slice or slice_counter is big enough to distribute the computation)
  // (so this doesn't finish the slice 
  while ((done_flag == false) && (slice_overflow_flag == false)) {
    
    
    // Loop through until we get to i=1 (so we defined a vector x)
    do {
      
      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }


      // Check to see if we should clear the slice_counter 
      // (i.e., when i=n since this means we're in a different slice)
      if (i == n) {
	cout << "slice_counter = " << slice_counter << endl;
	if (slice_counter >= SLICE_MIN)
	  slice_overflow_flag = true;

	slice_counter = 0;
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
    Q_val_double = precision() - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    Q_val = (unsigned long) round(Q_val_double);

    //    cout << " Float = " << Q_val_double << "   Long = " << Q_val << "  XX " << endl;
    /*
    cout << " The float value is " << Q_val_double << endl;
    cout << " The associated long value is " << Q_val << endl;
    cout << endl;
    */

    if (Q_val <= precision()) {
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
  }
	 
	 
  cout << "Leaving FastBinaryThetaInitial" << endl << endl;


  if (slice_overflow_flag == true)
    return x[n];   // Return the value of the slice we need to compute next!
  else
    return 1;      // same as above, but if we're done, then the value of the next slice should be 1.
                   // (even though we don't need to compute it!)
}




///////////////////////////////////////////////////////////////////////////////////////
// Reads the boolean ternary theta function for a slice from the temporary directory //
///////////////////////////////////////////////////////////////////////////////////////

bool boolean_ternary_theta::_read_from_distributed_host(const char* host) {
  //bool ReadTernaryThetaBinary_distributed_host(unsigned long theta_read[], long form[6], unsigned long precision, const char* host) {

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
	  DIST_DIR, host, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 


  // Check if the file exists
  if (FileExists(filename) == true) { 
    // Try to open the file, and return true if we're successful.
    cout << "Using filename " << filename << endl;
    cout << "Before read attempt" << endl;
    return read(filename);                 // This says when the read was successful.
    cout << "After read attempt" << endl;
  }
  else
    return false;  // This says we can't find the file. =(

 
  // Note: A false return value may mean that the file can't be found, or can't be read.


  // -------------------------------------------------------------------

  /*

  // Open the file for reading and check it opened correctly
  ifstream filein;
  filein.open(filename, ios::in | ios::binary);
  bool file_opened = filein.is_open();

  if (file_opened == true) {

    // Read the array
    filein.read( (char*) theta_read, sizeof(long) * precision());
    
    // Close the file
    filein.close();
  } 

  return file_opened;

  */

}




















// Main program (server) for computing the distributed theta function

void boolean_ternary_theta::compute_distributed(unsigned long SLICE_MIN) {
// void MakeTernaryThetaDistributed__by_x_n(long QQ[6], unsigned long Ternary_Precision) {

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

  cout << " Using SLICE_MIN = " << SLICE_MIN << endl;

  
  // Note: We distribute the computation after the slice size reaches SLICE_MIN


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
  

  // Start the computation for small slices
  long slice_todo;
  slice_todo = _compute_smallslices(Cholesky_temp, SLICE_MIN);
  //  slice_todo = FastBinaryThetaInitial(theta, Cholesky_temp, Ternary_Precision, SLICE_MIN);


  // DIAGNOSTIC:
  cout << " slice_todo = " << slice_todo << endl;
  cout << " The last two longs are: " << endl;
  cout << "  i = " << ((precision() >> 5) - 1) << " theta[i] = " << _theta[(precision() >> 5) - 1] << endl;
  cout << "  i = " << ((precision() >> 5) + 0) << " theta[i] = " << _theta[(precision() >> 5) + 0] << endl;
  cout << endl;



  // Check if we need to distribute
  if (slice_todo < 1) {

    cout << " Distributing the work..." << endl;

    // Allocate a temporary theta vector for the incoming slice information    
    boolean_ternary_theta theta_temp(QQ, precision());

        
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
	  hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(), slice_todo); 
      
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
    /*
    long i=0;
    while slice_to_do < 1
    */

    long i = 0;  
    bool all_slices_done_flag == false;
    vector<long> slices_to_redo;


    while ((num_of_hosts > 0) && (all_slices_done_flag == false)) {
      
      // Wait for a slice results file to appear, or a process to expire
      
      // INSERT: Sleep for 1 second.
      
      cout << "  Checking host #" << i << endl;
      
      // Check if a partial results file exists
      

      // If so, then incorporate the results
      if (theta_temp._read_from_distributed_host(hosts[i].c_str()) == true) {
      
	cout << " Recovered the data from host " << hosts[i] << endl;

	// Combine the new series with the old one
	for (unsigned long j = 0; j < precision(); j++)
	  _theta[j] = _theta[j] | theta_temp._theta[j];


	// Remove the temporary computation for that host
	char rm_filename[300];  
	sprintf(rm_filename, "rm -f %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
		DIST_DIR, hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	system("rm_filename");
	
	// Send another slice to that host             TODO: Check that there are still slices to send!!! =)
	char spawn_command[300];
	sprintf(spawn_command, "ssh %s 'nice ternary_theta_client %d %d %d %d %d %d %d %d'",
		hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(), slice_todo); 
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


	// Make the filename "Server_TEMPORARY__ternary_bool_theta__a_b_c_d_e_f__precision.bindata"
	char filename[200];  
	sprintf(filename, "%s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	
	// Write the temporary cumulatice theta function
	write(filename);


	// Make the filename "Server_TEMPORARY__ternary_bool_theta__a_b_c_d_e_f__precision.sliceinfo"
	//	char filename[200];  
	sprintf(filename, "%s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	

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
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(),
		DIST_DIR, "Server_CURRENT", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	system(mv_theta);
	char mv_sliceinfo[600];  
	sprintf(mv_sliceinfo, "mv -f %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(),
		DIST_DIR, "Server_CURRENT", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	system(mv_sliceinfo);

      }


      // Increment to cyclically loop through the hosts.
      if (i == num_of_hosts - 1)
	i = 0;
      else
	i++;

    }


    
    // Look for timeouts and re-delegate -- Still TO DO
    
    
    
  }
  else if (slice_todo > 1) 
    cout << " ERROR in MakeTernaryThetaDistributed__by_x_n: The current slice_todo " << slice_todo << " is > 1!" << endl;


  /*
  // Write the final output when we're done! =)
  cout << " about to write..." << endl;
  write();
  cout << " finished writing." << endl;
  */
}






// =========================================================================




// Main program (server) for computing the distributed theta function

void boolean_ternary_theta::compute_distributed_OLD(unsigned long SLICE_MIN) {
// void MakeTernaryThetaDistributed__by_x_n(long QQ[6], unsigned long Ternary_Precision) {

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

  cout << " Using SLICE_MIN = " << SLICE_MIN << endl;

  
  // Note: We distribute the computation after the slice size reaches SLICE_MIN


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
  

  // Start the computation for small slices
  long slice_todo;
  slice_todo = _compute_smallslices(Cholesky_temp, SLICE_MIN);
  //  slice_todo = FastBinaryThetaInitial(theta, Cholesky_temp, Ternary_Precision, SLICE_MIN);


  // DIAGNOSTIC:
  cout << " slice_todo = " << slice_todo << endl;
  cout << " The last two longs are: " << endl;
  cout << "  i = " << ((precision() >> 5) - 1) << " theta[i] = " << _theta[(precision() >> 5) - 1] << endl;
  cout << "  i = " << ((precision() >> 5) + 0) << " theta[i] = " << _theta[(precision() >> 5) + 0] << endl;
  cout << endl;



  // Check if we need to distribute
  if (slice_todo < 1) {

    cout << " Distributing the work..." << endl;

    // Allocate a temporary theta vector for the incoming slice information    
    boolean_ternary_theta theta_temp(QQ, precision());

        
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
	  hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(), slice_todo); 
      
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
    /*
    long i=0;
    while slice_to_do < 1
    */

    for(long i=0; i < num_of_hosts; i++) {

      cout << "  Checking host #" << i << endl;

      if (theta_temp._read_from_distributed_host(hosts[i].c_str()) == true) {

	cout << " Recovered the data from host " << hosts[i] << endl;

	// Combine the new series with the old one
	for (unsigned long j = 0; j < precision(); j++)
	  _theta[j] = _theta[j] | theta_temp._theta[j];


	// Remove the temporary computation for that host
	char rm_filename[300];  
	sprintf(rm_filename, "rm -f %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
		DIST_DIR, hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	system("rm_filename");
	
	// Send another slice to that host             TODO: Check that there are still slices to send!!! =)
	char spawn_command[300];
	sprintf(spawn_command, "ssh %s 'nice ternary_theta_client %d %d %d %d %d %d %d %d'",
		hosts[i].c_str(), QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(), slice_todo); 
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


	// Make the filename "Server_TEMPORARY__ternary_bool_theta__a_b_c_d_e_f__precision.bindata"
	char filename[200];  
	sprintf(filename, "%s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	
	// Write the temporary cumulatice theta function
	write(filename);


	// Make the filename "Server_TEMPORARY__ternary_bool_theta__a_b_c_d_e_f__precision.sliceinfo"
	//	char filename[200];  
	sprintf(filename, "%s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	

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
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(),
		DIST_DIR, "Server_CURRENT", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	system(mv_theta);
	char mv_sliceinfo[600];  
	sprintf(mv_sliceinfo, "mv -f %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo %s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.sliceinfo", 
		DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(),
		DIST_DIR, "Server_CURRENT", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
	system(mv_sliceinfo);

      }


      // Increment to cyclically loop through the hosts.
      if (i == num_of_hosts - 1)
	i = 0;
      else
	i++;

    }


    
    // Look for timeouts and re-delegate -- Still TO DO
    
    
    
  }
  else if (slice_todo > 1) 
    cout << " ERROR in MakeTernaryThetaDistributed__by_x_n: The current slice_todo " << slice_todo << " is > 1!" << endl;


  /*
  // Write the final output when we're done! =)
  cout << " about to write..." << endl;
  write();
  cout << " finished writing." << endl;
  */
}





/////////////////////////////////////////////////////////////////////////////////
// Test Code:
//
// 1. Check that MakeTernaryThetaDistributed__by_x_n agrees with 
//    MakeTernaryThetaDistributed when the slice_min is big, or the bound is small.
//



/*

bool TestDistribuedThetaComputation() {

  // Make the "Little Methuzalah" form
  // ans set the precision
  long QQ[6] = {2, 1, -1, 4, 3, 31};
  unsigned long precision = 100000;                    // *** Note: Ideally we'd be able to pass the precision, as well as the slice size!
  unsigned long theta_length = (precision >> 5) + 1;

  // Make two vectors to store the theta function
  unsigned long theta1[theta_length];
  unsigned long theta2[theta_length];


  // Make and read in the two forms
  MakeTernaryTheta(QQ,precision);
  ReadTernaryThetaBinary(theta1, QQ, precision);

  MakeTernaryThetaDistributed__by_x_n(QQ,precision);
  ReadTernaryThetaBinary(theta2, QQ, precision);

  
  // Compare them
  bool equal_flag = true;
  for(unsigned long i=0; i < theta_length; i++)
    if (theta1[i] != theta2[i]) {
      equal_flag = false;
      cout << "Error in TestDistribuedThetaComputation:" << endl;
      cout << "  The two little methuzalah theta functions differ in byte " << i << endl;
      cout << "    Undistributed: " << theta1[i] << "   Distributed: " << theta2[i] << endl;
      cout << "  Here are the previous entries: " << endl;
      cout << "    Undistributed: " << theta1[i-1] << "   Distributed: " << theta2[i-1] << endl;


    }
      
  return equal_flag;

}

*/








