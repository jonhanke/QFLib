


/////////////////////////////////////////////////////////
// This finds the ternary sublattice contained in the  //
// orthogonal complement of the upper left entry.      //
/////////////////////////////////////////////////////////

Matrix_mpz FindDecomposableTernarySublatticeComplementof1(Matrix_mpz Q) {

  // We assume that the upper-left entry is 1
  // and that Q is a 4 x 4.
  assert (Q(1,1) == 2);
  assert (Q.NumRows() == 4);
  assert (Q.NumCols() == 4);

  // For each row/column, perform elementary operations to cancel them out.
  for(long i=2; i<=4; i++){

    // Check if the (i,1)-entry is even (i.e. integral since we're uing 2Q here),
    // and double its row/ column if not.
    if (Q(i,1) % 2 != 0) 
      Q = MultiplySymmetric(Q, 2, i);

    // Now perform the (symmetric) elementary operations to cancel out the (i,1) entries/
    Q = AddSymmetric(Q, -Q(i,1)/2, i, 1);

  }

  // Check that we're done!
  if ((Q(1,2) != 0) || (Q(1,3) != 0) || (Q(1,4) != 0))
    cout << "ERROR IN FindDecomposableTernary.....of1 -- The form is not orthogonal... =(" << endl;
  
  // Return the lower left 3x3 block
  Matrix_mpz T;
  T.SetDims(3,3);
  for(long i=2; i<=4; i++)
    for(long j=2; j<=4; j++)
      T(i-1,j-1) = Q(i,j);

  return T;
}



//////////////////////////////////////////////////////////////////////
// Writes the theta function of a ternary form of desired precision //
//////////////////////////////////////////////////////////////////////

//! \deprecated  This has been moved to boolean_ternary_theta::compute().

void MakeTernaryTheta(long QQ[6], unsigned long Ternary_Precision) {

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


  // Use the Boolean_ternary_theta methods to do this
  boolean_ternary_theta temp_theta(QQ, Ternary_Precision);
  temp_theta.compute();
  temp_theta.write();


  // ---------------------------------------------------------------

  /*
  
  double Cholesky_temp[6];
  CholeskyTernary(QQ, Cholesky_temp);
  
  cout << " Computing the theta function " << endl;
  PrintTime();

  cout << "Allocating " << ((Ternary_Precision >> 5) +1) << " unsigned longs." << endl;
  //  unsigned long theta[(Ternary_Precision >> 5) +1];  // <-- This was the old allocation procedure... fails for 10^7.
  //  unsigned long *theta = (unsigned long *) calloc(sizeof(long), 10000000);
  unsigned long *theta = (unsigned long *) calloc(sizeof(long), ((Ternary_Precision >> 5) +1));
  cout << " Finished allocating the space!" << endl;

  FastBinaryTheta(theta, Cholesky_temp, Ternary_Precision);


  // DIAGNOSTIC:
  cout << " The last two longs are: " << endl;
  cout << "  i = " << ((Ternary_Precision >> 5) - 1) << " theta[i] = " << theta[(Ternary_Precision >> 5) - 1] << endl;
  cout << "  i = " << ((Ternary_Precision >> 5) + 0) << " theta[i] = " << theta[(Ternary_Precision >> 5) + 0] << endl;
  cout << endl;


  WriteTernaryThetaBinary(theta, QQ, Ternary_Precision); 

  */

}



//////////////////////////////////////////////////////////////////////
// Checks the representability of the square free numbers of a form //
//////////////////////////////////////////////////////////////////////

void Check4VarRepresentability(long Ternary[6], double B, long N, long char_top, 
			       unsigned long Precision, vector<long> & all_primes, 
			       Matrix_mpz Q, long form_num) {

  vector<long> Plist;
  vector<double> F4list;

  vector<mpz_class> Tlist;

  // Make a list of eligible primes and their values
  PrimeListFromBoundOptimized(Plist, F4list, B, N, char_top, all_primes); 
  cout << " There are " << Plist.size() << " eligible prime numbers." << endl;

  long VVnum = 5;
  if (Plist.size() < 5) 
    VVnum = Plist.size();
  if (Plist.size() > 0) {
    cout << " The last " << VVnum << " eligible primes and their F4(p) values are: " << endl;
    PrintTailV(Plist, VVnum);
    PrintTailV(F4list, VVnum);
  }
    


  // Read in the theta function
  boolean_ternary_theta bin_theta(Ternary, Precision);
  cout << "\n\n Started reading the ternary theta series." << endl;
  bin_theta.read();
  cout << " Finished reading the binary ternary theta series." << endl;




  // -----------------  Checking the squarefree numbers for the sublattice  ----------------------

  
  // Check the squarefree numbers
  vector<mpz_class> sublattice_exceptions;
  sublattice_exceptions = NewestSharpList(Tlist, B, Plist, F4list, bin_theta, 1, Q);


  // Write these numbers to a file

  char EXCEPTIONS_DIR[] = "/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/Exception_Data/";  // This is the global directory for the lists of exceptional numbers
  
  // Make the filename "Decomposable_sublattices/sublattice_of_form%d.txt"
  char filename[300];  
  sprintf(filename, "%s%ssublattice_of_form%d.txt", 
	  EXCEPTIONS_DIR, "Decomposable_sublattices/", form_num); 

  // Open the file for writing and check it opened correctly
  ofstream fileout;
  fileout.open(filename);

  if (! fileout.is_open())
    { cout << "Error opening output file " << filename << endl; exit (1); }

  // Write the list of sublattice exceptions
  for (unsigned long k=0; k < sublattice_exceptions.size(); k++)
    fileout << sublattice_exceptions[k] << endl;

  /*
  // Print the list of sublattice exceptions
  for (unsigned long k=0; k < sublattice_exceptions.size(); k++)
    cout << sublattice_exceptions[k] << ", ";
  */

  // Close the file
  fileout.close();

  // -------------------------------------------------------------------------------------------


  // Now check the exceptions of the sublattice ot find the exceptions of the original form.
  long Exception_Precision = sublattice_exceptions[sublattice_exceptions.size()-1].get_ui();  // Since it's already sorted.
  cout << "The largest sublattice exception is: " << Exception_Precision << endl;


  // ---------  Write the exceptions of the biglattice...  ---------------
    

  // Make the filename "Big_lattices/Form%d.txt"
  char filename1[300];  
  sprintf(filename1, "%s%sForm%d.txt", 
	  EXCEPTIONS_DIR, "Big_lattices/", form_num); 

  // Open the file for writing and check it opened correctly
  ofstream fileout1;
  fileout1.open(filename1);

  if (! fileout1.is_open())
    { cout << "Error opening output file " << filename1 << endl; exit (1); }

  
  // Find the excaptions of the big lattice if the sublattice candidates are < 10,000
  if (Exception_Precision < 10000){
    vector<bool> big_theta_series;

    // Make the filename
    char BIGLATTICE_THETA_DIR[] = "/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/Theta_Data/290_Biglattice_Thetas/";  // This is the global directory for the lists of exceptional numbers
    char filename[300];  
    sprintf(filename, "%sForm%d_Theta_10000.txt", BIGLATTICE_THETA_DIR, form_num);


    cout << "Reading in the theta series of the big form." << endl;
    
    // Read in the theta series of the big lattice
    ReadSeries_Boolean(big_theta_series, filename);

    cout << "Read in the theta series of the big form." << endl;

    // Check to see if the exceptions of the sublattice are represented by the lattice
    vector<mpz_class> biglattice_exceptions;
    for(unsigned long k=0; k < sublattice_exceptions.size(); k++) {
      mpz_class m;
      m = sublattice_exceptions[k];
      if (big_theta_series[m.get_ui()] == false) 
	biglattice_exceptions.push_back(m);
    }

    // Print the biglattice exceptions
    cout << "The Big lattice exceptions (from the sublattice) are:" << endl;
    for (unsigned long k=0; k < biglattice_exceptions.size(); k++)
      cout << biglattice_exceptions[k] << ", ";
    

    // Check something about their local conditions to make sure we don't miss any....



    // Write the list of sublattice exceptions
    for (unsigned long k=0; k < biglattice_exceptions.size(); k++)
      fileout1 << biglattice_exceptions[k] << endl;


  }
  else {
    cout << " Error: Some of the sublattice exceptions are bigger than 10,000. " << endl;
    fileout1 << "Error!  The candidates were too big!" << endl;
  }
  
  
  // Close the file
  fileout1.close();




}







//////////////////////////////////////////////////////////////////
// Identifies the upper-left ternary in each of our 6560 forms, //
// and describes those whose ternary is not regular.            //
//////////////////////////////////////////////////////////////////

void CheckAll6560Forms(vector<long> & all_primes) {

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
  /*  
  long JJ = 889;
  for(long k = JJ-1; k < JJ; k++) 
  */
  for(long k = 0; k < ListEnd; k++) 
    Theta_Precisions[k] = ThetaPrecision(bound_list[k], (lvl_list[k]).get_ui(), 
					 SquarefreePart(form_list[k].Determinant()).get_si(), 
					 1); // *** WARNING: This needs to be fixed!!! ***


  
  // Write the last 50 bounds
  cout << endl << endl;
  for(long k = ListEnd-50; k < ListEnd; k++) 
    cout << " The " << k+1 << "th bound is " << bound_list[sorted_indices[k]] 
	 << " which has Theta_precision " << Theta_Precisions[sorted_indices[k]] 
	 << endl;



  


  // Try the LocalSplitCovering Code on Form #223
  cout << "Finding the split covering for the form " << form_list[223] << endl;
  FindLocalSplitCovering(form_list[223]);
  cout << "Finished finding the split covering for the form " << form_list[223] << endl;





















  // Make the list of ternaries orthogonal to x^2.
  valarray <Matrix_mpz> orthogonal_ternary_list;
  orthogonal_ternary_list.resize(ListEnd);

  for(long k = 0; k < ListEnd; k++) 
    orthogonal_ternary_list[k] = FindDecomposableTernarySublatticeComplementof1(form_list[k]);


  /*
  while (1 == 1) {
  // Compute their theta functions to the appropriate precision
  cout << "\n\n Please select a group (0--65) of theta functions to compute: ";
  long i;
  cin >> i;
  
  long end;
  if (i+100 < ListEnd)
    end = 100*(i+1);
  else 
    end = ListEnd; 

  long tmp_form[6];
  Matrix_mpz tmp_T;
  for(long k = 100*i; k < end; k++) {
    tmp_T = orthogonal_ternary_list[k];
    tmp_form[0] = tmp_T(1,1).get_si() / 2;
    tmp_form[1] = tmp_T(1,2).get_si();
    tmp_form[2] = tmp_T(1,3).get_si();
    tmp_form[3] = tmp_T(2,2).get_si() / 2;
    tmp_form[4] = tmp_T(2,3).get_si();
    tmp_form[5] = tmp_T(3,3).get_si() / 2;

    MakeTernaryTheta(tmp_form, (unsigned long) ceil(Theta_Precisions[k]));
  }
  }
  */



  




  // Compute their theta functions to the appropriate precision
  cout << "\n\n Beginning to make ternary sublattice theta functions. " << endl;
  //  for(long k = 0; k < ListEnd; k++) {
  // /*
  for(long k = 6413; k < 6414; k++) {
  //for(long k = 0; k < 1000; k++) {
    long tmp_form[6];
    Matrix_mpz tmp_T;
    tmp_T = orthogonal_ternary_list[k];
    tmp_form[0] = tmp_T(1,1).get_si() / 2;
    tmp_form[1] = tmp_T(1,2).get_si();
    tmp_form[2] = tmp_T(1,3).get_si();
    tmp_form[3] = tmp_T(2,2).get_si() / 2;
    tmp_form[4] = tmp_T(2,3).get_si();
    tmp_form[5] = tmp_T(3,3).get_si() / 2;
    
    MakeTernaryTheta(tmp_form, (unsigned long) ceil(Theta_Precisions[k]));
  }
//  */
  cout << " Finished making ternary sublattice theta functions. \n\n" << endl;


  cout << "In the middle..." << endl;
  //  unsigned long ttt[10000000];
  cout << "In the middle..." << endl;


  // Start checking representability in order (((of ascending constants???)))  
  //  /*
  //  for(long k = 0; k < ListEnd; k++) {
  //  long KK = 75;
  long KK = 889;
  for(long k = KK-1; k < KK; k++) {
    long tmp_form[6];
    Matrix_mpz tmp_T;
    tmp_T = orthogonal_ternary_list[k];
    tmp_form[0] = tmp_T(1,1).get_si() / 2;
    tmp_form[1] = tmp_T(1,2).get_si();
    tmp_form[2] = tmp_T(1,3).get_si();
    tmp_form[3] = tmp_T(2,2).get_si() / 2;
    tmp_form[4] = tmp_T(2,3).get_si();
    tmp_form[5] = tmp_T(3,3).get_si() / 2;
    double Bound = bound_list[k];
    unsigned long lvl = lvl_list[k].get_ui();
    long char_top = SquarefreePart(form_list[k].Determinant()).get_si();
    unsigned long Precision = (unsigned long) ceil(Theta_Precisions[k]);

    cout << "\n\n\n Checking form #" << k+1 << endl;
    PrintTime();
    Check4VarRepresentability(tmp_form, Bound, lvl, char_top, Precision, all_primes, form_list[k], k+1);
  }
  //  */

}
