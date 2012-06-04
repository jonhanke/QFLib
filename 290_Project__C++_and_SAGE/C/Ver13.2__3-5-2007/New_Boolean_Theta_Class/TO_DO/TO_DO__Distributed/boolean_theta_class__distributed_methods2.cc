
///////////////////////////////////////////////////////////////////////
// Checks if a particular temporary computation results file exists
// (specified by a host and a slice)
///////////////////////////////////////////////////////////////////////


bool boolean_ternary_theta::_SliceFileExists(const long slice, const char* host) const {

  //  extern char DIST_DIR[];  // This is the global directory for the temporary distributed theta function files
    
  // Make the filename for each host "HOST__ternary_bool_theta__a_b_c_d_e_f__precision.bindata"
  char filename[300];  
  sprintf(filename, "%sSlice_%d__%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, slice, host, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 


  // Check if the file exists
  return FileExists(filename);
}




// ------------------------------------------------------------




///////////////////////////////////////////////////////////////////////////////////////
// Reads the boolean ternary theta function for a slice from the temporary directory //
///////////////////////////////////////////////////////////////////////////////////////

bool boolean_ternary_theta::_incorporate_temp_slice_results(const long slice, const char* host) {
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


  //  extern char DIST_DIR[];  // This is the global directory for the temporary distributed theta function files


  // Make the filename for each host "Slice_#__HOST__ternary_bool_theta__a_b_c_d_e_f__PRECISION.bindata"
  char filename[300];  
  sprintf(filename, "%sSlice_%d__%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, slice, host, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 

  
  // Allocate a temporary theta vector for the incoming slice information    
  boolean_ternary_theta theta_temp(QQ, precision());


  // Read the file in to the temporary theta vector
  // (TODO: Try to open the file, and return true if we're successful.)
  /*
  cout << "Using filename " << filename << endl;
  cout << "Before read attempt (in _incorporate_temp_slice_results())" << endl;
  */
  bool read_flag;
  read_flag = theta_temp.read(filename);     // NOTE: We should record this bit.  This says when the read was successful.
  /*
  cout << "After read attempt" << endl;
  */

  if (read_flag == true)
    cout << " Recovered the data from host " << host << endl;


  /*
  // DIAGNOSTIC -- see what we did
  for(long k = 0; k < precision(); k++)
    cout << " k = " << k << "    before = " << get_value(k) 
	 << "   new = " << theta_temp.get_value(k) 
	 << "   after = " << (get_value(k) | theta_temp.get_value(k)) << endl;
  */


  
  // Combine the new series with the old one
  for (unsigned long j = 0; j < _length(); j++)
    _theta[j] = _theta[j] | theta_temp._theta[j];




  // Remove the temporary computation we just read
  char rm_filename[310];  
  sprintf(rm_filename, "rm -f %s", filename); 
  system(rm_filename);
  

  return read_flag;  
}



/////////////////////////////////////////////////////
// Delegate a slice computation to a (remote) host //
/////////////////////////////////////////////////////

bool boolean_ternary_theta::_delegate_slice(distributed_info & D, long & slice, long i) const {

  //  extern char EXEC_DIR[];  // This is the global directory for the executable files
  /*  
  cout << " EXEC_DIR = " << EXEC_DIR << endl;
  */

  // Send another slice to that host             TODO: Check that there are still slices to send!!! =)
  char spawn_command[400];
  //  sprintf(spawn_command, "ssh -q -n %s 'nice %sternary_theta_client %d %d %d %d %d %d %d %d '  &",
  sprintf(spawn_command, "ssh -q -n %s 'nice %sternary_theta_client %d %d %d %d %d %d %d %d ' >> /dev/null &",
	  D.hosts[i].c_str(), EXEC_DIR, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision(), slice); 
  /*
  cout << "spawn_command = " << spawn_command << endl;
  */
  int return_code;
  return_code = system(spawn_command);


  if (true) {
    // If the delegation succeeds, then change that host's information...
    D.current_host_slice[i] = slice;
    time_t rawtime;
    time(&rawtime);
    D.host_assignment_time[i] = ctime(&rawtime);
    slice++;
    
    cout << " Sent slice " <<  D.current_host_slice[i] 
	 << " to " << D.hosts[i] << " at " << D.host_assignment_time[i] << endl; 
  }

  
  return true;  // TODO: Make this reflect the success of our delegation...
  
}





//////////////////////////////////////////////////////////////////////////
// Save the current state in two temporary files (theta and state file) //
//////////////////////////////////////////////////////////////////////////

bool boolean_ternary_theta::_save_current_distributed_state(const distributed_info & D) const {

  //  extern char DIST_DIR[];  // This is the global directory for the temporary distributed theta function files
    
  // Make the filename "Server_TEMPORARY__ternary_bool_theta__a_b_c_d_e_f__precision.bindata"
  char filename[200];  
  sprintf(filename, "%s%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, "Server_TEMPORARY", QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
  
  // Write the temporary cumulative theta function
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
  for (long j=0; j < D.num_of_hosts ; j++) 
    fileout << D.hosts[j] << "  /  " 
	    << D.current_host_slice[j] << "  /  " 
	    << D.host_assignment_time[j] << "  /  " << endl;

  
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


  /*
  double Cholesky_temp[6];
  CholeskyTernary(QQ, Cholesky_temp);
  */



  // Start the computation for small slices
  long slice_todo;
  slice_todo = _compute_smallslices(SLICE_MIN);
  //  slice_todo = FastBinaryThetaInitial(theta, Cholesky_temp, Ternary_Precision, SLICE_MIN);


  /*
  // DIAGNOSTIC -- see what we did in compute_smallslices()
  cout << "\n\n Results from the small slices..." << endl;
  for(long k = 0; k < precision(); k++)
    cout << " k = " << k << "    value = " << get_value(k) << endl;
  */  

  /*
  // DIAGNOSTIC:
  cout << " slice_todo = " << slice_todo << endl;
  cout << " The last two longs are: " << endl;
  cout << "  i = " << ((precision() >> 5) - 1) << " theta[i] = " << _theta[(precision() >> 5) - 1] << endl;
  cout << "  i = " << ((precision() >> 5) + 0) << " theta[i] = " << _theta[(precision() >> 5) + 0] << endl;
  cout << endl;
  */


  // Check if we need to distribute
  if (slice_todo < 1) {

    cout << " Distributing the work..." << endl;
        
    // Make a list of the available computers
    vector<string> hosts_list;

    // These are the test machines
    hosts_list.push_back("tux1");
    hosts_list.push_back("tux2");
    hosts_list.push_back("tux3");


    /*
      // These are the Duke Machines
      hosts_list.push_back("tux1");
      hosts_list.push_back("tux2");
      hosts_list.push_back("tux3");
      hosts_list.push_back("tux4");
      hosts_list.push_back("tux5");
      hosts_list.push_back("tux6");
      hosts_list.push_back("tux7");
      hosts_list.push_back("tux8");
      hosts_list.push_back("tux9");
    */

    // Instantiating (making) the distributed info class with our hosts. =)
    distributed_info Distributed(hosts_list);

    cout << " There are " << Distributed.num_of_hosts << " hosts." << endl;

    //    char hosts[num_of_hosts][100] = 

        // Some other possibilities:
        // -------------------------
        //     vector<char[100]> hosts(num_of_hosts)= {"tux1", "tux2", "tux3"};
        //
        // vector<string> hosts(num_of_hosts);
        // hosts = {"tux1", "tux2", "tux3"};



    // Make a directory to look for the partial results
    extern char DIST_DIR[];  // This is the global directory for the temporary distributed theta function files


    // Delegate chunks in order
    for(long i=0; i < Distributed.num_of_hosts; i++) {
      
      if (_delegate_slice(Distributed, slice_todo, i) == false) 
	cout << " Error in Delegating slice " << slice_todo << endl;
      
      
    }
    
    
    
    // Look for results and re-delegate
    /*
      long i=0;
      while slice_to_do < 1
    */
    
    long i = 0;  
    bool all_slices_done_flag = false;
    vector<long> slices_to_redo;
    

    while ((Distributed.num_of_hosts > 0) && (all_slices_done_flag == false)) {
      
      // Wait for a slice results file to appear, or a process to expire
      
      // /*
      // INSERT: Sleep for 1 second.
      system("sleep 1; date");
      cout << "  Checking host #" << i << endl;
      // */
      
      
      // Check if a partial results file exists
      if (_SliceFileExists(Distributed.current_host_slice[i], Distributed.hosts[i].c_str()) == true) {
	
	// If so, try to incorporate the results.  
	_incorporate_temp_slice_results(Distributed.current_host_slice[i], Distributed.hosts[i].c_str());

	// If there are more slices to delegate, try to delegate them
	if (slice_todo != 1) {

	  
	  if (_delegate_slice(Distributed, slice_todo, i) == false) 
	    cout << " Error in Delegating slice " << slice_todo << endl;
	  

	}
	else {
	  Distributed.current_host_slice[i] = 1;
	  // If there are no more slices, then decrease the number of active calculations.
	}

	// Finally, save the state of the computation. (one of every 100 times)
	if ((slice_todo % 100 == 0) || (slice_todo == 1))
	  _save_current_distributed_state(Distributed);
		
      }
      else {
	// Either there is no partial results file, (so it's still computing)
	// or there was a failure in reading it...
	//
	// If this fails, then add this slice to the failure list and remove the host...
      }

      
      
      // Increment to cyclically loop through the hosts.
      if (i == Distributed.num_of_hosts - 1)
	i = 0;
      else
	i++;


      // Test to see if we're done with all computations
      all_slices_done_flag = true;
      for(long j=0; j< Distributed.num_of_hosts; j++)
	if (Distributed.current_host_slice[j] != 1)
	  all_slices_done_flag = false;  
      
    }


    
    // Look for timeouts and re-delegate -- Still TO DO
    
    
    
  }
  else if (slice_todo > 1) 
    cout << " ERROR in MakeTernaryThetaDistributed__by_x_n: The current slice_todo " << slice_todo << " is > 1!" << endl;


  //  /*
  // Write the final output when we're done! =)
  cout << " about to write..." << endl;
  write();
  cout << " finished writing." << endl;
  // */
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








