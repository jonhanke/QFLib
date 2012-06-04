

//////////////////////////////////////////
// Read in the list of 34 ternary forms //
//////////////////////////////////////////

void Read34Ternaries(long QQ[34][6]) {
  long form_count = 0;
  char c;
  
  // Make the filename
  char ternaryfilename[100];  
  sprintf(ternaryfilename, "/home/postdoc/jonhanke/290_Project/Ternaries.txt"); 
  
  // Open the file for writing and check it opened correctly
  ifstream ternaryfile;
  ternaryfile.open(ternaryfilename);
  
  if (! ternaryfile.is_open())
    { cout << "Error opening output file " << ternaryfilename; exit(1); }
  
  
  // Main loop -- repeat until we find a leading "]"
  do {
    
    // Read until we find "["
    do {
      ternaryfile >> c;
    } while (c != '[');
    
    // Read in the 9 matrix entries
    mpq_class x[9];
    for (long i = 0; i<9; i++) {    
      ternaryfile >> x[i];  // Read the entry
      ternaryfile >> c;     // Read "," or "]"
    }  
    
    // Print the input
    /*
      cout << endl;
      cout << "c = " << c << endl;
      cout << "[ " << x[0] << " " << x[1] << " " << x[2] << " ]" << endl;
      cout << "[ " << x[3] << " " << x[4] << " " << x[5] << " ]" << endl;
      cout << "[ " << x[6] << " " << x[7] << " " << x[8] << " ]" << endl;
      cout << endl;
    */
    
    // Test the input
    bool ok_flag = true;
    if (c != ']') 
      ok_flag = false;  // Check that there are only 6 entries.
    if ((x[1] != x[3]) || (x[2] != x[6]) || (x[5] != x[7]))
      ok_flag = false;  // Check the matrix is symmetric
    if ((x[0].get_den() != 1) || (x[4].get_den() != 1) || (x[8].get_den() != 1))
      ok_flag = false;  // Check the diagonal is integral
    
    
    // Double the strictly upper triangular entries
    x[1] = 2*x[1];
    x[2] = 2*x[2];
    x[5] = 2*x[5];
    x[1].canonicalize();
    x[2].canonicalize();
    x[5].canonicalize();
    
    
    // Test the strictly upper triangular entries are integral
    if ((x[1].get_den() != 1) || (x[2].get_den() != 1) || (x[5].get_den() != 1))
      ok_flag = false;  
        
    // Abort if the tests fail
    if (ok_flag == false) {
      cout << " There is a problem with a form..." << endl;
      exit(1);
    } 
    /*
      else 
      cout << " Successfully read form " << form_count << "... =)" << endl;
    */
    
    
    // Add the form to the list 
    QQ[form_count][0] = x[0].get_num().get_si();
    QQ[form_count][1] = x[1].get_num().get_si();
    QQ[form_count][2] = x[2].get_num().get_si();
    QQ[form_count][3] = x[4].get_num().get_si();
    QQ[form_count][4] = x[5].get_num().get_si();
    QQ[form_count][5] = x[8].get_num().get_si();
    
    // Increment the form counter
    form_count++;
    
    // Read the next chacter to see if it's a ","
    ternaryfile >> c;   
    
  } while (c == ',');
  
  
    
  // Print the 34 forms we just read
  cout << "\nHere are the 34 forms we read: " << endl;
  for (long i = 0; i < 34; i++) {
    cout << "  " << i+1 << ":  [ ";
    for (long j = 0; j < 6; j++)
      cout << QQ[i][j] << "  ";
    cout << "]" << endl;
  }
}  




////////////////////////////////////////////////////////////////////////
// Makes any of the 34 ternary theta functions with precision <= 10^7 //
////////////////////////////////////////////////////////////////////////

void MakeOneThetaFrom34Ternaries() {

  // Read in the 34 ternary forms
  long QQ[34][6];
  Read34Ternaries(QQ);
  
  // Find the level of each of the 34 ternary forms
  for(long i=0; i<34; i++)
    cout << " The level of form #" << i+1 << " is " << QF_Ternary_Level(QQ[i]) << endl;
  
  // Choose a theta function...
  cout << "Enter a form (0--33) to compute the theta function of: ";
  long i;
  cin >> i;


  // Do the computation ans write the result
  const unsigned long Ternary_Precision = 10000000;  
  boolean_ternary_theta theta3(QQ[i], Ternary_Precision);
  
  cout << " Computing Form #" << i << endl;
  PrintTime();

  theta3.compute();
  theta3.write();
}




// Look for the exceptions < 10^7 for the 34 ternaries
void FindTernaryExceptions_by_theta() {
    
  // Read in the 34 ternary forms
  long QQ[34][6];
  Read34Ternaries(QQ);
  
  // Find the level of each of the 34 ternary forms
  for(long i=0; i<34; i++)
    cout << " The level of form #" << i+1 << " is " << QF_Ternary_Level(QQ[i]) << endl;
  
  
  // Look for the exceptions < 10^7 for the 34 ternaries
  const unsigned long Ternary_Precision = 10000000;
  for (long i = 19; i < 34; i++) {
    cout << " Checking Form #" << i << endl;
    cout << " [ " << QQ[i][0] << ", "
	 << QQ[i][1] << ", " << QQ[i][2] << ", "
	 << QQ[i][3] << ", " << QQ[i][4] << ", "
	 << QQ[i][5] << " ] " << endl;
    PrintTime();

    
    // Read in the ternary theta series
    boolean_ternary_theta theta3(QQ[i], Ternary_Precision);
    cout << "Started reading the Theta series" << endl;
    theta3.read();
    cout << "Finished reading the Theta series" << endl;

    
    // Check the local conditions
    long local_mod;
    vector<long> missed_classes;
    vector<long> local_mod_vector;
    
    
    cout << "Finding the local conditions" << endl; 
    
    FindTernaryLocalConditions(local_mod, missed_classes, QQ[i], local_mod_vector); 
    
    cout << "Found the local conditions" << endl; 
    
    /*
      cout << "local_mod = " << local_mod << endl;
      cout << "local_mod / 4 = " << (local_mod/4) << endl;
      cout << "missed_classes = ";
      for (unsigned long j = 0; j < missed_classes.size(); j++)
      cout << missed_classes[j] << ", ";
    */  
    
    // Modify the local_mod_vector at p=2 to exclude multiples of mod/4 or mod for each prime power.
    local_mod_vector[0] = local_mod_vector[0] / 4; 
    
    
    // Look for the exceptions
    cout << " The exceptions for form " << i << " are: " << endl;
    cout << " local_mod_vector.size() = " << local_mod_vector.size() << endl;
    cout << " We exclude exceptional numbers divisible by: " << endl;
    for (unsigned long k=0; k < local_mod_vector.size(); k++)
      cout << local_mod_vector[k] << ", ";
    cout << endl << endl;
    
    for(unsigned long j=0; j < Ternary_Precision; j++) {
      if (theta3.get_value(j) == false) {
	
	// Check if it is locally missed
      	bool is_missed = false;
	for (unsigned long k=0; k < missed_classes.size(); k++)
	  if ((j % local_mod) == (unsigned long) missed_classes[k])
	    is_missed = true;
	
	// Eventually we'll check if it is divisible by the square of
	// an anisotropic prime, and then divide it out and check again...
	// BUT FOR NOW JUST CHECK IT'S RELATIVELY PRIME TO local_mod...
	// -------------------------------------------------------------
	// ACTUALLY, WE'LL JUST EXCLUDE THE NUMBERS DIVISIBLE BY each prime power local_mod/4
	// and local_mod/4 when p=2.
	for (unsigned long k=0; k < local_mod_vector.size(); k++)
	  if ((j % local_mod_vector[k]) == 0) 
	    is_missed = true;
	
	
	// Write the exceptions
	if (is_missed == false)
	  cout << " " << j << endl;
      }
    }
    cout << endl << endl;
  }
}  





// Look for the exceptions < 10^7 for the 34 ternaries
void FindTernaryExceptions_New() {
    
  // Read in the 34 ternary forms
  long QQ[34][6];
  Read34Ternaries(QQ);
  
  // Find the level of each of the 34 ternary forms
  for(long i=0; i<34; i++)
    cout << " The level of form #" << i+1 << " is " << QF_Ternary_Level(QQ[i]) << endl;
  
  
  // Look for the exceptions < 10^7 for the 34 ternaries
  const unsigned long Ternary_Precision = 10000000;
  for (long i = 0; i < 34; i++) {
    cout << " Checking Form #" << i+1 << endl;
    cout << " [ " << QQ[i][0] << ", "
	 << QQ[i][1] << ", " << QQ[i][2] << ", "
	 << QQ[i][3] << ", " << QQ[i][4] << ", "
	 << QQ[i][5] << " ] " << endl;
    PrintTime();

    
    // Read in the ternary theta series
    boolean_ternary_theta theta3(QQ[i], Ternary_Precision);
    cout << "Started reading the Theta series" << endl;
    theta3.read();
    cout << "Finished reading the Theta series" << endl;

    
    // Check the local conditions
    vector<long> local_mod_vector;
    long local_repn_array[200][9];  // NOTE: We should really not have a 200 here.  This should be dynamic, and depend on the number of p|N.
    vector<long> aniso_vector;
    
    cout << "Finding the local conditions" << endl; 
    
    FindTernaryLocalConditions_New(QQ[i], local_mod_vector, local_repn_array, aniso_vector); 
    
    cout << "Found the local conditions" << endl; 
    
    /*
      
      cout 
      cout << "local_mod = " << local_mod << endl;
      cout << "local_mod / 4 = " << (local_mod/4) << endl;
      cout << "missed_classes = ";
      for (unsigned long j = 0; j < missed_classes.size(); j++)
      cout << missed_classes[j] << ", ";
    */  
    
    // Modify the local_mod_vector at p=2 to exclude multiples of mod/4 or mod for each prime power.
    local_mod_vector[0] = local_mod_vector[0] / 4; 
    
    
    // Look for the exceptions
    cout << " The exceptions for form " << i << " are: " << endl;
    cout << " local_mod_vector.size() = " << local_mod_vector.size() << endl;
    cout << " These are the local prime power moduli used in checking local representability: " << endl;
    for (unsigned long k=0; k < local_mod_vector.size(); k++)
      cout << local_mod_vector[k] << ", ";
    cout << " These are the local representability conditions (array): " << endl;
    for (unsigned long k=0; k < local_mod_vector.size(); k++) {
      for (unsigned long l=0; l < 9; l++)
	cout << local_repn_array[k][l] << ", ";
      cout << endl;
    }
    cout << " There are " << aniso_vector.size() << " anisotropic primes. " << endl;
    cout << " The anisotropic primes are: " << endl;
    for (unsigned long k=0; k < aniso_vector.size(); k++) {
	cout << aniso_vector[k] << ", ";
    }
    cout << endl << endl;

    for(unsigned long j=0; j < Ternary_Precision; j++) {
      if (theta3.get_value(j) == false) {


	//	cout << "\n Checking " << j << endl;

	// Divide out all even powers of anisotropic primes
	long m = j;
	for (unsigned long k=0; k < aniso_vector.size(); k++)
	  if (m % (aniso_vector[k] * aniso_vector[k]) == 0)
	    m = m / (aniso_vector[k] * aniso_vector[k]);

	//	cout << " after removing anisotropic squares, we get " << m << endl;


	// Check for local representability at each eligible prime
	bool is_locally_missed = false;
	for (unsigned long k=0; k < local_mod_vector.size(); k++) {  // NOTE: Could speed this up by exiting as soon as is_missed is true
	
	  long p = local_repn_array[k][0];
	  unsigned long v = Valuation(m, p);

	  //	  cout << " Hi -- Using m = " << m << " and p = " << p << endl;
	  //	  cout << "  valuation = " << v << endl;
	  
	  
	  // Find the associated number kk (so m = x * p^(2*kk))
	  unsigned long kk = v / 2;
	  long m1 = m / LongPow(p, 2*kk);

	  //	  cout << "  kk = " << kk << endl;  


	  // Find the squareclass index (with 1 <= index <= 8  or  1 <= index <= 4)
	  long index = 0;
	  if (p > 2) {
	    if (m1 % p == 0) {
	      m1 = m1 / p;
	      index = 2;
	    }

	    if (LegendreSymbol(mpz_class(m1), mpz_class (p)) == 1)
	      index = index + 1;
	    else 
	      index = index + 2;
	  }	  

	  // ( Treat p = 2 separately )
	  else {                       
	    if (m1 % p == 0) {
	      m1 = m1 / p;
	      index = 4;
	    }

	    if (m1 % 8 == 1)
	      index = index + 1;
	    else if (m1 % 8 == 3)
	      index = index + 2;
	    else if (m1 % 8 == 5)
	      index = index + 3;
	    else if (m1 % 8 == 7)
	      index = index + 4;
	  }

	  //	  cout << "  index = " << index << endl;

	  
	  // Check if j (or m) is locally represented
	  if ( (local_repn_array[k][index] < 0) || ((long) kk < local_repn_array[k][index])) {
	    is_locally_missed = true;
	    //	    cout << "   not locally repn at p = " << p << endl;
	  }
	  //	  else 
	  //	    cout << "   locally repn at p = " << p << endl;
	}
	
	
	// Write the exceptions
	if (is_locally_missed == false)
	  cout << " " << j << endl;
      }
    }
    cout << endl << endl;
  }
}  







//////////////////////////////////////////////////////////////////
// Identifies the upper-left ternary in each of our 6560 forms, //
// and describes those whose ternary is not regular.            //
//////////////////////////////////////////////////////////////////

void Find_4var_with_Irreg_Ternaries() {


  // List the forms which don't have any of the regular ternaries inside,
  // and find (and sort) their associated constants.
  long regular_ternary_list[24] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 22, 24, 28, 30, 34};
  //long irregular_ternary_list[10] = {13, 18, 23, 25, 26, 27, 29, 31, 32, 33};
  vector<long> irregular_quaternary_list; 

  // Read in the 34 ternary forms
  long Ternary[34][6];
  Read34Ternaries(Ternary);

  
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
  cout << "\n Reading the 290 cusp data file " << endl;
  Read290Data("/home/postdoc/jonhanke/290_Project/290_Cusp_info/290-cusp-all.m", form_list, lvl_list, cusp_const_list);
  cout << "\n Finished reading the 290 cusp data file \n" << endl;


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
  /*
  cout << "They are: " << endl;
  for (unsigned long k = 0; k < unproven_form_number_list.size(); k++)
    cout << unproven_form_number_list[k] << ", ";

  cout << "\n They have bounds: " << endl;
  for (unsigned long k = 0; k < unproven_form_number_list.size(); k++)
    cout << bound_list[unproven_form_number_list[k]] << ", ";
  */
  cout << "\n The big bounds are: " << endl;
  for (unsigned long k = 0; k < unproven_form_number_list.size(); k++)
    if (bound_list[unproven_form_number_list[k]] > 500)
      cout << bound_list[unproven_form_number_list[k]] << ", ";


  //  /*

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

  //  */




  // Make the ThetaPrecision vectors
  //  ThetaPrecision(double B, long N, long chi_top, long diag_coeff)
  vector<double> Theta_Precisions;
  Theta_Precisions.resize(ListEnd);    
  for(long k = 0; k < ListEnd; k++) 
    Theta_Precisions[k] = ThetaPrecision(bound_list[k], (lvl_list[k]).get_ui(), 
					 SquarefreePart(form_list[k].Determinant()).get_si(), 
					 1); // *** WARNING: This needs to be fixed!!! ***


  
  // Write the bounds
  cout << endl << endl;
  for(long k = ListEnd-50; k < ListEnd; k++) 
    cout << " The " << k+1 << "th bound is " << bound_list[sorted_indices[k]] 
	 << " which has Theta_precision " << Theta_Precisions[sorted_indices[k]] 
	 << endl;



  /*
  sort(sorted_bounds.begin(), sorted_bounds.end());

  cout << endl << endl;
  for(long k = ListEnd-10; k < ListEnd; k++) 
    cout << " The " << k+1 << "th bound is " << sorted_bounds[k] << endl;
  */





  /*
    void sort_example() 
    // illustrates the use of the sort algorithm
    // see alg7.cpp for complete source code
    {
    // fill both a vector and a deque
    // with random integers
    std::vector<int> aVec(15);
    std::deque<int> aDec(15);
    std::generate(aVec.begin(), aVec.end(), randomValue);
    std::generate(aDec.begin(), aDec.end(), randomValue);
    
    // sort the vector ascending
    std::sort(aVec.begin(), aVec.end());
    
    // sort the deque descending
    std::sort(aDec.begin(), aDec.end(), greater<int>() );
    
    // alternative way to sort descending
    std::sort(aVec.rbegin(), aVec.rend());
    }    
  */

}


