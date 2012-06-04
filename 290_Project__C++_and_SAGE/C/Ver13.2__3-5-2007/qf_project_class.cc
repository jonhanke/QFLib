
// Constructor to make a project directory system (if it doesn't exist) 
QF_Project::QF_Project(const string & projectname, const string & filename, 
		       const string & cusp_const_dir, const string & cusp_const_prefix, const set<long> & excluded_forms) 
  :QF_Datafiles(projectname) {


  // Make the list of excluded forms
  Excluded_forms = excluded_forms;

  
  /*
  // Initialize the data folders and names for this project 
  QF_Datafiles::QF_Datafiles(projectname);
  */


  // Make the form list filename  --  "~/QF_Project_Data/projectname/projectname__form_list.txt"
  string absolute_form_list_filename, absolute_cusp_const_dir;
  absolute_form_list_filename =  Project_Dir + projectname + "__form_list.txt";
  absolute_cusp_const_dir =  Project_Dir + "Cuspidal_Constants/";


  // Make the form list file if there isn't one in the project directory 
  if (FileExists(absolute_form_list_filename) == false) {
    
    if ((filename == "") || (FileExists(GetAbsolutePath(filename)) == false)) {
      string command_line;
      command_line = "echo '[ ]' >  " + absolute_form_list_filename;
      system(command_line.c_str());
    } else  {                        
      string command_line;
      command_line = "cp " + GetAbsolutePath(filename) + " " + absolute_form_list_filename;
      system(command_line.c_str());
    }   

  }


  // Make the Cuspidal Constant Directory if there isn't one in the project directory 
  if (FileExists(absolute_cusp_const_dir) == false) {
    
    if ((cusp_const_dir == "") || (FileExists(GetAbsolutePath(cusp_const_dir)) == false)) {
      cout << "Error in QF_Project Cosntructor:  Cusp_Const_Dir doesn't make sense!  Using the directory " << GetAbsolutePath(cusp_const_dir) << endl;
      assert(0==1);
    } else  {                        
      string command_line;
      command_line = "cp -rp " + GetAbsolutePath(cusp_const_dir) + " " + absolute_cusp_const_dir;
      system(command_line.c_str());
    }   

  }


  // Read in the forms from the datafile  -- FIX ALL OF THIS!!!
  /*
  Read_Approx_Forms_and_Cusp_Bounds(absolute_form_list_filename);
  */
  Read_Exact_Forms_and_Cusp_Bounds(absolute_form_list_filename, absolute_cusp_const_dir, cusp_const_prefix);


  // Compute their levels and store them too
  for (long i=0; i<Form_List.size(); i++) {
    //cout << i << "  " << Form_List[i];        // DIAGNOSTIC
    Level_List.push_back(Form_List[i].QFLevel());
  }

  // Allocate space for the exception lists.
  Exception_Lists.resize(Form_List.size());

}




/////////////////////////////////////////////////////////////////////////////////////////////////
// Reads William Stein's Magma formatted list of forms and approximate constants -- DEPRECATED //
/////////////////////////////////////////////////////////////////////////////////////////////////

void QF_Project::Read_Approx_Forms_and_Cusp_Bounds(const string & filename) {

  // Try opening data file
  ifstream datafile;
  datafile.open(filename.c_str(), ios::in);
  if (! datafile.is_open())
    { cout << "Error opening file"; exit (1); }

  
  // Declare some variables
  char c;
  int num;
  mpz_class lvl;
  mpq_class temp_entry;
  Matrix_mpz M(4,4);
  float cusp_const, error_const;


  // Begin reading data:
  // -------------------

  // Eat the initial "["
  datafile >> c;
  //    cout << " 1 - Read in " << c << endl;


  // Read in all of the data... =)   <------  NEED TO FIX THIS TO DETERMINE THE NUMBER OF FORMS FIRST, BY PRE-READING!
  for(size_t k=0; k<6560; k++) {

    /*
    // Skip some lines
    cout << "\n\n";
    cout << " Starting with k = " << k << endl;
    */

    // Eat the initial "<"
    datafile >> c;
    //      cout << " 2 - Read in " << c << endl;


    // ---------------------------------------------

    // Read in the Matrix           <------  EVENTUALLY WANT TO FIX THIS TO TAKE ANY SIZE MATRICES...
    for(size_t i=1; i<=4; i++)
      for(size_t j=1; j<=4; j++) {
	
	// Eat the "[" or ","
	datafile >> c;
	//          cout << " 3 - Read in " << c << endl;

	// Read in the number of the quadratic form
	datafile >> temp_entry;
	M(i,j) = 2*temp_entry;
	
      }

    // Add the matrix to the list
    Form_List.push_back(M);

    /*
    cout << "The matrix is \n" << endl;
    cout << M << endl;
    */
    
    // Eat the closing "],"
    datafile >> c;
    //      cout << " 4 - Read in " << c << endl;
    datafile >> c;
    //      cout << " 5 - Read in " << c << endl;


    // ---------------------------------------------

    // Read the cusp constant
    datafile >> cusp_const;
    Cusp_Constant_List.push_back(cusp_const);
    //    cout << " The cusp constant is " << cusp_const << endl;


    // Eat the ">," or ">]"
    datafile >> c;
      //cout << " Read in " << c << endl;
    datafile >> c;
      //cout << " Read in " << c << endl;
    
    
  }

  
}




///////////////////////////////////////////////////////////////
// Reads the list of forms and the exact cuspidal constants. //
///////////////////////////////////////////////////////////////

void QF_Project::Read_Exact_Forms_and_Cusp_Bounds(const string & filename, const string & cusp_dir, const string & cusp_prefix) {

  // Read in the list of forms
  //string new_filename = "104_auxiliary_quaternaries.txt";
  cout << "Using filename = " << filename << endl;
  string new_filename = filename;


  // Try opening the data file with the list of forms (as on the 290 website)
  ifstream datafile;
  datafile.open(new_filename.c_str(), ios::in);
  if (! datafile.is_open())
    { cout << "Error opening file"; exit (1); }

  
  // Declare some variables
  char c;
  int num;
  mpz_class lvl;
  mpq_class temp_entry;
  Matrix_mpz M(4,4);
  float cusp_const, error_const;
  
  
  // Begin reading data:
  // -------------------
  
  // Eat the initial "["
  datafile >> c;
  //cout << " 1 - Read in " << c << endl;
  
  
  
  bool list_end_flag = false;
  size_t k = 1;
  while (list_end_flag == false) {
    
    // Read in the Matrix           <------  EVENTUALLY WANT TO FIX THIS TO TAKE ANY SIZE MATRICES...
    for(size_t i=1; i<=4; i++)
      for(size_t j=1; j<=4; j++) {
	
	// Eat the "[" or ","
	datafile >> c;
	//cout << " 2 - Read in " << c << endl;
	
	// SANITY CHECK: Check that we haven't ended prematurely, and the format is ok. =)
	assert((c == '[') || (c == ','));
	
	// Read in the coefficient of the quadratic form
	datafile >> temp_entry;
	//cout << " Read temp_entry as " << temp_entry << endl;
	M(i,j) = 2*temp_entry;
	
      }
    
    // Add the matrix to the list
    Form_List.push_back(M);


    // DIAGNOSTIC
    /*
    cout << "The matrix is \n" << endl;
    cout << M << endl;
    */

    
    // Eat the closing "],"
    datafile >> c;
    //cout << " 3 - Read in " << c << endl;
    datafile >> c;
    //cout << " 4 - Read in " << c << endl;
    
    
    // SANITY CHECK: Check that we haven't ended prematurely, and the format is ok. =)
    assert((c == ']') || (c == ','));
    
    
    // Detect the end of the list
    if (c != ',') 
      list_end_flag = true;
    else
      k++;
    
    
    // DIAGNOSTIC
    if (list_end_flag == true)
      cout << " Read in " << k << " quadratic forms from " << new_filename <<"! =)" << endl;

    
    
    // ---------------------------------------------
    
    /*
    // Read the cusp constant
    datafile >> cusp_const;
    Cusp_Constant_List.push_back(cusp_const);
    //    cout << " The cusp constant is " << cusp_const << endl;
    
    
    // Eat the ">," or ">]"
    datafile >> c;
    //cout << " Read in " << c << endl;
    datafile >> c;
    //cout << " Read in " << c << endl;
    */
    
  }
  
  
  
  // Populate the cuspidal constant list with -1's by default:
  // ---------------------------------------------------------
  cout << endl << " Starting to populate the cusp constant list with -1's" << endl;
  for(size_t k=0; k < Form_List.size(); k++) 
    Cusp_Constant_List.push_back(-1);
  cout << " Finished populating the cusp constant list with -1's" << endl << endl;



  // Run through the forms looking for a cuspidal constant file for each one! =) 
  // ---------------------------------------------------------------------------
  /*
  string cuspidal_constant_folder = "Auxiliary_Constants/";
  string cusp_file_prefix = "Aux_const_";
  */
  string cuspidal_constant_folder = cusp_dir;
  string cusp_file_prefix = cusp_prefix;
  string cusp_file_suffix = ".txt";
  char tmp_char;
  float tmp_cusp_const;

  cout << "Starting to read the cuspidal constant summary files." << endl;

  // Read in all of the cusp form constants (and check that the forms agree with our list)
  for(size_t k=0; k < Form_List.size(); k++) {

    // DIAGNOSTIC
    /*
    cout << endl << " Starting to read the cusp constant for form #" << k+1 << endl;
    */

    // Make the filename
    char num_data[10];
    sprintf(num_data, "%lu", k+1);
    string new_cusp_filename = cuspidal_constant_folder + cusp_file_prefix + string(num_data) + cusp_file_suffix;

    
    // Try opening the cusp constant file for form number k, and if we fail, assign it's constant as -2 (to indicate the error).
    ifstream datafile;
    datafile.open(new_cusp_filename.c_str(), ios::in);
    if (! datafile.is_open()) { 
      // cout << "Error opening file " << new_cusp_filename << " for form #" << k+1 << endl;   
	Cusp_Constant_List[k] = -2; 
    }
    else {
        
      // Try to read the cuspidal constant
      datafile >> tmp_cusp_const;
      //cout << " Read " << tmp_cusp_const;
      
      // Read the rest of the file and compare that with the form listed in Form_list

      /*
      // Read the opening '['
      datafile >> tmp_char;
      cout << " Read " << tmp_char;
      assert(tmp_char == '[');
      */
      
      // Read in the matrix of the form corresponding to this constant
      for(size_t i=1; i<=4; i++)
	for(size_t j=1; j<=4; j++) {
	  
	  // Eat the "[" or ","
	  datafile >> tmp_char;
	  //cout << " 2 - Read in " << tmp_char << endl;
	  
	  // SANITY CHECK: Check that we haven't ended prematurely, and the format is ok. =)
	  assert((tmp_char == '[') || (tmp_char == ','));
	  
	  // Read in the coefficient of the (already doubled) quadratic form
	  datafile >> temp_entry;
	  //cout << " Read temp_entry as " << temp_entry << endl;
	  M(i,j) = temp_entry;
	  
	}
      
      // Read the closing ']'
      datafile >> tmp_char;
      //cout << " Read " << tmp_char;
      assert(tmp_char == ']');
      
      
      // SANITY CHECK: Make sure they agree
      /*
	cout << " M = " << M << endl;
	cout << " Form_List[k] = " << Form_List[k] << endl;
      */
      assert(M == Form_List[k]);
      
      // If they agree, then add the constant to the list.
      Cusp_Constant_List[k] = tmp_cusp_const;

    }
      
  }


  // Status report
  cout << "Finished reading the cuspidal constant summary files." << endl;
  
}




//////////////////////////////////////////////////////////////////
// Identifies the upper-left ternary in each of our 6560 forms, //
// and describes those whose ternary is not regular.            //
//////////////////////////////////////////////////////////////////

void QF_Project::FindExceptions(const string & type, const long & first, const long & last) {

  /*
  // List the forms which don't have any of the regular ternaries inside,
  // and find (and sort) their associated constants.
  long regular_ternary_list[24] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 22, 24, 28, 30, 34};
  //long irregular_ternary_list[10] = {13, 18, 23, 25, 26, 27, 29, 31, 32, 33};
  vector<long> irregular_quaternary_list; 

  // Read in the 34 ternary forms
  long Ternary[34][6];
  Read34Ternaries(Ternary);
  */

  /*  
  // Prepare to read in the 6560 forms
  valarray <Matrix_mpz> form_list;
  valarray <mpz_class> lvl_list;
  valarray <float> cusp_const_list;

  valarray <mpq_class> local_const_list;
  valarray <float> computed_lower_bound_list;
  valarray <mpq_class> approximate_lower_bound_list;

  const long ListEnd = 6560;  // This is the number of quaternary forms to read.

  form_list.resize(ListEnd);
  lvl_list.resize(ListEnd);
  cusp_const_list.resize(ListEnd);


  // Read in the 6560 quadratic forms
  cout << "\n Reading the 290 cusp data file ---------------- "; 
  PrintTime();
  Read290Data("/home/postdoc/jonhanke/290_Project/290_Cusp_info/290-cusp-all.m", form_list, lvl_list, cusp_const_list);
  cout << "\n Finished reading the 290 cusp data file ------- ";
  PrintTime();
  cout << endl << endl;


  // Read in the 6560 F4 Constants
  valarray <double> bound_list;  
  bound_list.resize(ListEnd);
  const char BOUND_DIR[] = "/home/postdoc/jonhanke/290_Project/F4_Constants/";

  for(long k = 0; k < ListEnd; k++) {
    
    // Make the filename "/home/postdoc/jonhanke/290_Project/F4_Constants/
    //                    FormXXXX_F4_bound.txt"
    char filename[200];  
    sprintf(filename, "%sForm%d_F4_bound.txt", BOUND_DIR, k+1); 
    
    
    // Open the file for reading and check it opened correctly
    ifstream filein;
    filein.open(filename);
    
    if (! filein.is_open())
      { cout << "Error opening input file " << filename; exit (1); }
    
    // Read the array
    filein >> bound_list[k];
    
    // Close the file
    filein.close();
  }
  */
 
  /*
  // Extract the upper left ternary

  long extracted_ternary_list[ListEnd][6];
  for(long k = 0; k < ListEnd; k++) {
    extracted_ternary_list[k][0] = form_list[k](1,1).get_si() / 2;
    extracted_ternary_list[k][1] = form_list[k](1,2).get_si();
    extracted_ternary_list[k][2] = form_list[k](1,3).get_si();
    extracted_ternary_list[k][3] = form_list[k](2,2).get_si() / 2;
    extracted_ternary_list[k][4] = form_list[k](2,3).get_si();
    extracted_ternary_list[k][5] = form_list[k](3,3).get_si() / 2;
  }
  cout << " Finished extracting upper left ternaries. " << endl;


  // Check each form for a regular ternary, and list those without one.
  cout << " Begin checking the extracted ternaries. " << endl;
  vector<long> unproven_form_number_list;
  valarray <double> extracted_ternary_num_list;  
  extracted_ternary_num_list.resize(ListEnd);
  for(long k = 0; k < ListEnd; k++) {
    //         cout << "  Starting k = " << k << endl;
    bool done_flag = false; 
    long ternary_index = 0;

    // Loop through all 34 ternaries looking for a match
    while ((done_flag == false) && (ternary_index < 34)) {
      bool match_flag = true;
      for(long l = 0; l < 6; l++) 
	if (extracted_ternary_list[k][l] != Ternary[ternary_index][l])
	  match_flag = false;

      ternary_index ++;

      if (match_flag == true) {
	done_flag = true;
	ternary_index --;
      }

    }

    // If we found one of the 34 ternaries
    if (done_flag == true) {

      // Check if the ternary is regular
      bool regular_flag = false;
      for(long l = 0; l < 24; l++)
	if (ternary_index == regular_ternary_list[l])
	  regular_flag = true;
      
      if (regular_flag == false)
	unproven_form_number_list.push_back(k);      // Note: The numbering starts form zero here...
    }
    
    // if we didn't find a matching ternary then give an error message 
    else
      cout << "ERROR in identifying ternaries:  Ternary not found... =(" << endl;
  }
  cout << " Finished checking the extracted ternaries. " << endl;


  // Run through the unproven form list
  cout << endl << endl;
  cout << "There are " << unproven_form_number_list.size() << " forms (with irregular ternaries) to consider." << endl;
  cout << "\n The big bounds are: " << endl;
  for (unsigned long k = 0; k < unproven_form_number_list.size(); k++)
    if (bound_list[unproven_form_number_list[k]] > 500)
      cout << bound_list[unproven_form_number_list[k]] << ", ";
*/



  /*

  // Sort the bounds by their constants
  vector<double> sorted_bounds;
  vector<long> sorted_indices;      // Staring at zero.

  for(long k = 0; k < ListEnd; k++) 
    sorted_bounds.push_back(bound_list[k]);
  for(long k = 0; k < ListEnd; k++) 
    sorted_indices.push_back(k);

  bool swap_flag = true;
  while (swap_flag == true) {
    swap_flag = false;
    
    // Swap adjacent unordered entries
    for(long k = 0; k < ListEnd-1; k++) {
      if (sorted_bounds[k] > sorted_bounds[k+1]) {
	double tmp_bound;
	tmp_bound = sorted_bounds[k];
	sorted_bounds[k] = sorted_bounds[k+1];
	sorted_bounds[k+1] = tmp_bound;
	
        long tmp_index;
        tmp_index = sorted_indices[k];
        sorted_indices[k] = sorted_indices[k+1];
        sorted_indices[k+1] = tmp_index;
	
	swap_flag = true;
      }
    }
  }

  */



  /*

  // Make the ThetaPrecision vectors
  //  ThetaPrecision(double B, long N, long chi_top, long diag_coeff)
  vector<double> Theta_Precisions;
  Theta_Precisions.resize(ListEnd);    
  for(long k = 0; k < ListEnd; k++) 
    Theta_Precisions[k] = ThetaPrecision(bound_list[k], (lvl_list[k]).get_ui(), 
					 SquarefreePart(form_list[k].Determinant()).get_si(), 
					 1); // *** WARNING: This needs to be fixed!!! ***

  */


  
  /*

  // Write the last few F4 bounds
  long last_nn = 10;
  cout << endl << endl;
  cout << " Printing the largest " << last_nn << " F4 bounds: " << endl;
  cout << " ------------------------------------------" << endl;
  for(long k = max(0, ListEnd - last_nn); k < ListEnd; k++) 
    cout << " The " << k+1 << "th bound is " << bound_list[sorted_indices[k]] 
	 << " which has Theta_precision... ??? " 
	 << endl;
  cout << endl << endl;

  */




  /*

  // Testing with the 420 form:
  // --------------------------
  Matrix_mpz Q420(4,4);
  Q420(1,1) = 2;
  Q420(2,2) = 6;
  Q420(3,3) = 10;
  Q420(4,4) = 14;


  // Check the first 100 local densities
  Q420.CheckSiegelRange(1, 100);
  

  // Finding the Eisenstein coefficient lower bound
  cout << " Attempting to find the Eisenstein coefficient lower bound:" << endl;
  mpq_class E420_bound, E420_num_bound;
 
  E420_bound = Q420.EisensteinLowerBound();
  cout << " The theoretical lower bound is: " << E420_bound << endl << endl;
  
  E420_num_bound = Q420.NumericalEisensteinLowerBound();
  cout << " The numerical lower bound is: " << E420_num_bound << endl << endl; 


  // Summary
  cout << " The theoretical lower bound is: " << E420_bound << " which is about " << (E420_bound.get_d() * 1.0) << endl;
  cout << " The numerical lower bound is: " << E420_num_bound << " which is about " << (E420_num_bound.get_d() * 1.0) << endl;
  cout << endl;

  // Sanity check
  if (E420_bound <= E420_num_bound)
    cout << " The theoretical lower bound is less than the numerical lower bound from the coefficients! =)" << endl;
  else
    cout << " ERROR: There's a problem since the theoretical lower bound is larger than what we observe in the coefficients! =( " << endl;

  cout << " =============================================================================================" << endl << endl;
  

  */


  /*
  // Testing with form #224:
  // -----------------------
    
  Matrix_mpz Q;
  long form_num = 224;  // This is the number of the form (starting from #1 and ending at #6560)
  Q = form_list[form_num - 1];
  cout << " Checking Form #" << form_num << ": " << endl;
  cout << " ==================== " << endl << endl;
  
    
  // Check the first 100 local densities
  Q.CheckSiegelRange(1, 100);  
    
    
  // Determine the exceptions of Q
  Representability Q_Exceptions(Q, cusp_const_list[form_num - 1]);
  cout << " The (square-free) exceptions of Q are: " <<  Q_Exceptions.GetExceptions() << endl;
  */


  // --------------------------------------------------------------------------------------------


  cout << " There are " << Form_List.size() << " forms to check the exceptions of..." << endl << endl;


  vector<long> overflow_list;    // List of forms for which Representability gives and overflow


  // Run through all forms:
  // ----------------------
  //for(long k=0; k<Form_List.size(); k++) {

  //  for(long k=0; k<6560; k++) {
  // for(long k=0; k<100; k++) {
  // for(long k=46; k<47; k++) {              // Diagnostic for Form #47
  // for(long k=250; k<251; k++) {              // Diagnostic for Form #251
  // for(long k=3494; k<3495; k++) {              // Diagnostic for Form #3495
  // for(long k=4309; k<4310; k++) {              // Diagnostic for Form #3410
  // for(long k=6229; k<6230; k++) {              // Diagnostic for Form #6230


  //  Distributed Computation Loops:
  //  ------------------------------
  //  for(long k=0; k<500; k++) {             // This is for Tux1
  //  for(long k=500; k<1000; k++) {          // This is for Tux2
  //  for(long k=1000; k<1500; k++) {         // This is for Tux3
  //  for(long k=1500; k<2000; k++) {         // This is for Tux4
  //    for(long k=2000; k<2500; k++) {         // This is for Tux5
  //  for(long k=2500; k<3000; k++) {         // This is for Tux6
  //  for(long k=3000; k<3500; k++) {         // This is for Tux7
  //  for(long k=3500; k<4000; k++) {         // This is for Tux8
  //  for(long k=4000; k<4500; k++) {         // This is for Tux9
  //  for(long k=4500; k<5000; k++) {         // This is for Navier
  //  for(long k=5000; k<5500; k++) {         // This is for Stokes
  //  for(long k=5500; k<6000; k++) {         // This is for Austin
  //  for(long k=6000; k<6500; k++) {         // This is for Navier2
  //    for(long k=6500; k<6560; k++) {         // This is for Austin2

  for(long k = first-1; k <= last-1; k++) {      // Run through the specified range of forms 


    // Only compute forms not explicitly excluded! =)
    if (Excluded_forms.find(k+1) == Excluded_forms.end()) {     
      
      Matrix_mpz Q;
      Q = Form_List[k];
      cout << endl << endl;
      cout << " Checking Form #" << k+1 << ": " << endl;
      cout << " ==================== " << endl << endl;
      PrintTime();    
      cout << endl;
      
      
      /*
      // Check the theta series (Magma vs. Computed)
      Q.Check_ComputeTheta_vs_MagmaTheta();
      
      // Check the first 100 local densities
      Q.CheckSiegelRange(1, 100, (*this).Eis_Dir);  
      */
      
      
      // Determine the exceptions of Q
      Representability Q_Exceptions(*this, Q, Cusp_Constant_List[k], type);
      if (Q_Exceptions.overflow == false) {
	cout << " The (square-free) exceptions of Q are: " << Q_Exceptions.GetExceptions() << endl;
	cout << " Adding these to the exception list! " << endl;
	vector<mpz_class> tmp_new_exceptions;
	tmp_new_exceptions = Q_Exceptions.GetExceptions();
	Project_exception_set.insert(tmp_new_exceptions.begin(), tmp_new_exceptions.end());
	Exception_Lists[k] = tmp_new_exceptions;   // Add these exceptions to the overall project list
      }
      else {
	cout << " *** OVERFLOW FOR THIS FORM! *** " << endl;
	cout << " Adding it to the overflow list! " << endl;
	overflow_list.push_back(k);
	assert(0==1);
      }
      
      
      // Print a running total of the (squarefree) exceptions
      cout << " The running list of squarefree exceptions is: " << endl << Project_exception_set << endl;
      
    }
    else 
      cout << "\n\n\n *************** EXCLUDED FORM # " << k+1 << " ********************\n\n\n";
    
    
  }
  

  cout << " ------------------------------------------------------------------------------------------- " << endl ;  
  cout << " ------------------------------------------------------------------------------------------- " << endl << endl << endl;  


  // Print the overflow forms
  cout << " The forms (#1-" << Form_List.size() << ") which overflowed are: " << endl << "   ";
  if (overflow_list.size() >= 1)
    for (long k=0; k < overflow_list.size(); k++)
      cout << (overflow_list[k] + 1) << ", ";
  else
    cout << "No Overflow Forms! =)" << endl;
  cout << endl << endl;


  /*
  // Sort the exception list
  if (exception_list.size() >= 2) {

    bool swap_flag2 = true;
    while (swap_flag2 == true) {
      swap_flag2 = false;
      
      // Swap adjacent unordered entries
      for(long k = 0; k < exception_list.size() - 1; k++) {
	if (exception_list[k] > exception_list[k+1]) {
	  mpz_class tmp;
	  tmp = exception_list[k];
	  exception_list[k] = exception_list[k+1];
	  exception_list[k+1] = tmp;
	  
	  swap_flag2 = true;
	}
      }
    }
    
  }
  */


  cout << " =============================================================================================" << endl;
  cout << "                                Finished Computing the Exceptions!" << endl;
  cout << " =============================================================================================" << endl << endl << endl;


  
}




// Print the exceptions for the current range of forms
void QF_Project::PrintExceptions(const long & first, const long & last) {

  // Local Variables
  long number_of_exception_forms = 0;

  // Say which forms have exceptions
  cout << " Non-universal (non-excluded) forms from Form #" << first << "--" << last << ":" << endl; 
  cout << " -------------------------------------------------------------- " << endl << endl;

  // Loop through all of the forms in our specified range
  for(long k = first-1; k <= last-1; k++) {    

    // Only compute forms not explicitly excluded! =)
    if (Excluded_forms.find(k+1) == Excluded_forms.end()) 

      // Check if there are any exceptions
      if (Exception_Lists[k].empty() == false) {
	cout << " Form #" << (k+1) << " has exceptions " << Exception_Lists[k] << endl;
        number_of_exception_forms++;

      }
  }


  // Print the project exception list
  cout << " The exceptions are: " << endl << "   ";
  if (Project_exception_set.size() >= 1) 
    for (set<mpz_class>::iterator iter=Project_exception_set.begin(); iter != Project_exception_set.end(); iter++)
      cout << *iter << ", ";
  else
    cout << "No Exceptions! =)" << endl;
  cout << endl << endl << endl;


  // Print the excluded form list
  cout << " The excluded forms are: " << endl << "   ";
  if (Excluded_forms.size() >= 1) 
    for (set<long>::iterator iter=Excluded_forms.begin(); iter != Excluded_forms.end(); iter++)
      cout << "#" << *iter << ", ";
  else
    cout << "No Excluded Forms! =)" << endl;
  cout << endl << endl << endl;


  // Print the final statistics
  cout << " Looking from Form #" << first << "--" << last << endl << endl; 
  cout << " There were " << Excluded_forms.size() << " excluded forms." << endl;
  cout << " There were " << number_of_exception_forms << " forms which were not universal." << endl;
  cout << " There were a total of " << Project_exception_set.size() << " exceptions for the project." << endl;


}
