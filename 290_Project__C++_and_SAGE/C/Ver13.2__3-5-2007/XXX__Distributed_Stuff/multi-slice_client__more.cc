


////////////////////////////////////////////////////////////////
// Writes the boolean ternary theta function as a binary file //
// NOTE: This is just a copy of the routine in maketheta.cc!  //
// (with some minor modifications to the directory and filename...)
////////////////////////////////////////////////////////////////

void boolean_ternary_theta::_WriteTernaryThetaBinary_client_new(char* hostname) {

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


  // Make the filename for each host "Done__HOST__ternary_bool_theta__a_b_c_d_e_f__PRECISION.bindata"
  char filename[300];  
  sprintf(filename, "%sDone__%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, hostname, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 


  // Make the filename for each host "TEMP__Done__HOST__ternary_bool_theta__a_b_c_d_e_f__PRECISION.bindata"
  char temp_filename[300];  
  sprintf(temp_filename, "%sTEMP__Done__%s__ternary_bool_theta__%d_%d_%d_%d_%d_%d__%u.bindata", 
	  DIST_DIR, hostname, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 


  // Write the temporary file
  write(temp_filename);


  // Move the temporary file to it's real location
  char mv_command[600];  
  sprintf(mv_command, "mv %s %s", temp_filename, filename); 
  system(mv_command);

}







// Read the vector of slices for this host
vector<long> boolean_ternary_theta::_read_slicevector(const char* host) {

  //  extern char DIST_DIR[];  // This is the global directory for the temporary files

  vector<long> slice_vec;
  char c;
  long l;


  // Make the filename for the slice vector 
  char filename[200];  
  sprintf(filename, "%s%s__%d_%d_%d_%d_%d_%d__%u.slicevectordata", 
	  DIST_DIR, host, QQ[0], QQ[1], QQ[2], QQ[3], QQ[4], QQ[5], precision()); 
  
  // Open the vector for writing
  ifstream filein;
  filein.open(filename);
  
  if (! filein.is_open())
    { cout << "Error opening output file " << filename; exit (1); }


  // Read the slice vector to a file for that host
  filein >> c;                                  // Reads the '['
  filein >> c;                                  // Checks if we have a number (by pre-reading a char)
  while (c != ']') {
    filein.putback(c);        // Put back the pre-read character
    filein >> l;
    slice_vec.push_back(l);
    filein >> c;               // Check out the next character
  }
  

  // Close the file
  filein.close();


  // DIAGNOSTIC
  cout << "slice_vec.size() = " << slice_vec.size() << endl;
  cout << " Read the vector of slices: " << endl;
  cout << "[ ";
  for (long j=0; j < slice_vec.size(); j++) 
    cout << slice_vec[j] << " " << endl;
  cout << "]";


  // Return the vector of slices
  return slice_vec;
}








// Main program (server) for computing the distributed theta function

void boolean_ternary_theta::compute_client_new() {

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


  // Read in the slice vector
  vector<long> slice_vec;
  slice_vec = _read_slicevector(hostname);


  // Start the computation for these slice
  cout << "\n\nStarting the slice" << endl << endl;
  _compute_with_slices(Cholesky_temp, 0, slice_vec);
  cout << "\n\nFinishing the slice" << endl << endl;



  // DIAGNOSTIC -- to check the computed values before writing
  for(long kk=0; kk < precision(); kk++)
    cout << " kk =  " << kk << "    value = " << get_value(kk) << endl;


  // Write the output to the Distributed directory
  _WriteTernaryThetaBinary_client_new(hostname); 

}








