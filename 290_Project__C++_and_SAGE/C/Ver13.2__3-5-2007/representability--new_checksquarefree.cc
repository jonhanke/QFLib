
////////////////////////////////////////////////////////////////////////////////////////////
// Incrementally generates eligible squarefree numbers, and checks their representability.
/////////////////////////////////////////////////////////////////////////////////////////////

void Representability::_CheckSquarefree_New(const QF_Datafiles & qf_data) {


  // Make a counter to count squarefree numbers
  // and lists of exceptions and overflows...
  mpz_class Counter;
  Counter = 0;
  //  vector<mpz_class> exception_list, overflow_list, overflow_differences, overflow_depths;



  int Max_num = _max_num_prime_factors;

  vector<mpz_class> T_list;


  // Stores the possible exceptions for each # of prime factors
  vector< vector<mpz_class> > possible_exception_list_vec, possible_exception_depth_vec;    
  possible_exception_list_vec.resize(_max_num_prime_factors+1);    // We index this with the number of prime factors (which is >=1).
  possible_exception_depth_vec.resize(_max_num_prime_factors+1);   // We index this with the number of prime factors (which is >=1).


  long cover_number = 1;  // Keeps track of the local cover we're using to check the forms.
  


  // Check if the big form is locally universal (once and for all)
  bool Bigform_universal_flag = _local_cover.Big_Local().IsUniversal();




  // Get its boolean ternary theta functions (relative to the maximum allowed depth)
  cout << "   There are " << _local_cover.multiplicity()  << " split local coverings." << endl;
  cout << "   Computing the relevant boolean ternary theta series for local cover #" << cover_number << ": " << endl;
  _Get_Cover_Ternary_Thetas(qf_data, cover_number);
  cout << "   Finished computing the relevant boolean ternary theta series for local cover #" << cover_number << "! " << endl << endl << endl;
  




  

  // Find the local conditions for QQ

  // TO DO:
  // ------
  // This needs to be a valarray<mpz_class>...
  //  vector<long> aniso_vector = _local_cover.BigForm.AnisotropicPrimes();
  vector<long> aniso_vector;

  /*
  // Some old notation:
  // ------------------

  QQ_sub = QQ as a split form:

  QQ =  [ d  0 ]        d = number
        [ 0  T ]        T = ternary form

  Ternary = T 

  theta3 = boolean_ternary_theta
  */




    /*

    // Print the local conditions
    cout << "Here are the local conditions: " << endl;
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

    */
  


  // ======================================================================================================


  /*
  // DIAGNOSTIC
  cout << " The list of eligible primes is " << endl << _plist << endl;
  */


  // ======================================================================================================

  //////////////////////////////////////////////////////////////
  // Make the array of Max_num pointers and step through them //
  // to get a minimal list of square-free numbers t to check. //
  //////////////////////////////////////////////////////////////

    //  long PP_length = _plist.size();
  int carry_ptr;
  int counter = 0;

  cout << "\n Starting to check number with " << _max_num_prime_factors << " prime factors. "  << endl;
  cout << "   ";
  PrintTime();


  // Loop through all the possible #'s of prime factors
  for (int num = _max_num_prime_factors; num >= 1; num--) {          
    
    // Initialize the array for j prime factors
    int Last_ptr = num - 1;
    carry_ptr = Last_ptr;  // This keeps the compiler from complaining... =)
    vector<long> Ptr_array(num);
    for (int i=0; i<= Last_ptr; i++) 
      Ptr_array[i] = i;

    // Step through the possible prime factors and construct eligible squarefree numbers
    bool done_flag = false;
    while (done_flag == false) {

      bool carry_flag = false;    // Start without wanting to carry.

      // Check that we're still in the allowed primes range
      if (Ptr_array[Last_ptr] < _plist.size()) { 
            
	// Compute F4(t) for the product of the currently indexed primes
	double temp_F4_prod = 1;
	for (int i=0; i<=Last_ptr; i++)
	  temp_F4_prod = temp_F4_prod * _f4list[Ptr_array[i]];  

	// If F4(t) < bound and it's locally represented, add t to the list and reset the carry_ptr 
	if (temp_F4_prod <= _F4_bound) {
	  mpz_class temp_squarefree = 1;
	  for (int i=0; i<=Last_ptr; i++) 
	    temp_squarefree = temp_squarefree * _plist[Ptr_array[i]];
	  
	  /*
	  // DIAGNOSTIC
	  if (temp_squarefree == 290)
	     cout << " *************** Adding 290 to the list! ************************* " << endl;
	  */

	  // DIAGNOSTIC
	  //	  cout << "About to check the local representability of " << temp_squarefree << endl;
	  
	  // Pre-Check that the number is locally represented (or BigForm is locally universal)
	  // (NOTE: It is probably faster to check the local representability *after* it is found as an exception...) 
	  if ((Bigform_universal_flag == true) || 
	      (_local_cover.Big_Local().IsLocallyRepresented__Squarefree(temp_squarefree) == true)) {

	    // DIAGNOSTIC
	    //	    cout << temp_squarefree << " is locally represented! =)" << endl;

	    T_list.push_back(temp_squarefree);
	    //        T_file << temp_squarefree << ", ";
	    Ptr_array[Last_ptr] = Ptr_array[Last_ptr] + 1;
	    carry_ptr = Last_ptr;
	  }
	  else 
	    carry_flag = true;
	  
	} 
	else	
	  carry_flag = true;
      } 
      else 
	carry_flag = true;
      
      // If we're out of range or F4(t) is too big, then try to carry
      if (carry_flag == true) {

	// Perform the carry if we're not at the first array entry
	if (carry_ptr > 0){
	  carry_ptr = carry_ptr - 1;
	  Ptr_array[carry_ptr] = Ptr_array[carry_ptr] + 1;
	  for (int i=carry_ptr + 1; i<=Last_ptr; i++)
	    Ptr_array[i] = Ptr_array[carry_ptr] + (i - carry_ptr);
	}

	// If we're at the first entry, then we're done. =)
	else
	  done_flag = true;
      }


      // Check the numbers in Tlist if there are enough, or if we're done.
      if ((T_list.size() == 1000000) || (done_flag == true)) {
	Counter = Counter + T_list.size();

	/*
	// OVERFLOW DIAGNOSTIC...
	T_list.resize(0);
	T_list.push_back(mpz_class(940869930));
	*/

	/*
	// DIAGNOSTIC
	cout << "Made the list of square-free numbers: " << endl << T_list << endl;
	*/



	cout << "      About to check " << Counter << " numbers " << endl << endl;




	/*
	// Check these squarefree numbers for execptions
	_NewestCheckNumbers_by_cover(qf_data, T_list, CHECKING_DEPTH);	
	*/




	// Check the vector of numbers, and add them to the current list
	vector<mpz_class> new_possible_exception_list;
	vector<mpz_class> new_possible_exception_depths;
	long num_of_prime_factors = Last_ptr+1;      // This is passed to tell the local cover statistics class how to count these numbers!
    	_NewestCheckNumbers_by_cover_Subroutine(qf_data, T_list, num_of_prime_factors, new_possible_exception_list, new_possible_exception_depths);

	possible_exception_list_vec[num_of_prime_factors].insert(possible_exception_list_vec[num_of_prime_factors].end(), 
								 new_possible_exception_list.begin(), new_possible_exception_list.end());
	possible_exception_depth_vec[num_of_prime_factors].insert(possible_exception_depth_vec[num_of_prime_factors].end(),
								  new_possible_exception_depths.begin(), new_possible_exception_depths.end());





	T_list.resize(0);

	cout << "    Cumulative results: " << endl;
	cout << "      Checked " << Counter << " numbers using the first local cover" << endl;
	if (new_possible_exception_list.size() > 0) {
	  cout << "      There are " << new_possible_exception_list.size() << " new possible exceptions." << endl;
	  cout << "      They are: " << new_possible_exception_list << endl;
	}
	cout << "      Number of possible exceptions for " << num_of_prime_factors << " prime divisors: " 
	     << possible_exception_list_vec[num_of_prime_factors].size() << endl;
	//	cout << "      Number of exceptions: " << _squarefree_exceptions.size() << endl;
	//	cout << "      Number of overflows: " << _squarefree_overflows.size() << endl;
	cout << "      "; PrintTime();
	cout << endl;
      }
	

    }


    // /*
    // Description of the number of square-free numbers with num prime factors:
    //    cout << "\n Decreasing the number of prime factors." << endl;
    //    cout << "\n   Last_ptr being decreased from " << Last_ptr << " to " << Last_ptr - 1 << endl;
    //    cout << "   Added " << T_list.size() - counter << " square-free numbers with " 
    //   	 << num << " prime factors." << endl;
    //    cout << "   The T_list has " << T_list.size() << " numbers." << endl;
    cout << "\n Decreasing the number of prime factors from " << Last_ptr + 1 << " to " << Last_ptr << endl;
    cout << "   ";
    PrintTime();

    counter = T_list.size();  // Reset the counter for the next round. 
    // */


  }





  // Set the largest exception allowed for checking with the Quaternary_Exceptions class
  const long BIGGEST_EXCEPTION_FOR_EASY_CHECKING = 10000;
  mpz_class max_eligible = 0;

  // Sort the exception lists (for each number of prime factors)
  vector<mpz_class> tmp_max_eligible_vec;
  tmp_max_eligible_vec.resize(_max_num_prime_factors);
  for(long ind=1; ind < _max_num_prime_factors; ind++) {
    sort(possible_exception_list_vec[ind].begin(), possible_exception_list_vec[ind].end());
    if (possible_exception_list_vec[ind].empty() == false)
      tmp_max_eligible_vec[ind] = possible_exception_list_vec[ind][possible_exception_list_vec[ind].size() - 1];
  }

  // Find the largest eligible exception (if there are some)            
  for(long ind=1; ind < _max_num_prime_factors; ind++)
    max_eligible = max(max_eligible, tmp_max_eligible_vec[ind]);


  // If any of the remaining possible exceptions are too big, and we don't have too many local covers, try another local cover!
  if (max_eligible > BIGGEST_EXCEPTION_FOR_EASY_CHECKING) {


    // Status report after that checking...
    cout << endl;
    cout << endl;
    cout << "The statistics for checking with this local cover (#" << cover_number << ") are:" << endl;
    cout << endl;
    cout << _local_stats_vec[cover_number-1] << endl;    // Note: We start indexing from 0 here.
    cout << endl;
    for(long ind=1; ind <= _max_num_prime_factors; ind++) {
      cout << " There are " << possible_exception_list_vec[ind].size() << " remaining possible exceptions." << endl;
      cout << "   They are: " << possible_exception_list_vec[ind] << endl;
      cout << endl;
    }
    cout << " The largest possible exception " << max_eligible << " is larger than " << BIGGEST_EXCEPTION_FOR_EASY_CHECKING << endl;
    cout << " so we're going to try using another local cover!" << endl;
    cout << endl;


    // DIAGNOSTIC -- EXIT!
    cout << " ERROR! We are trying to use only use one local cover!" << endl;
    exit(1);


    // Make another local cover
    cover_number++;
    if ( cover_number > _local_cover.multiplicity())
      _local_cover.IncrementLocalSplitCoverings();

    // Get its boolean ternary theta functions (relative to the maximum allowed depth)
    cout << "   There are " << _local_cover.multiplicity()  << " split local coverings." << endl;
    cout << "   Computing the relevant boolean ternary theta series for local cover #" << cover_number << ": " << endl;
    _Get_Cover_Ternary_Thetas(qf_data, cover_number);
    cout << "   Finished computing the relevant boolean ternary theta series for local cover #" << cover_number << "! " << endl << endl << endl;
    

    /*   ****** FIX THIS! ******
    // Use this to check the remaining numbers
    T_list = possible_exception_list_vec[XXXXX];
    possible_exception_list.clear();
    long num_of_prime_factors = 0;   // This says that we store the remaining numbers as a big bunch, regardless of the number of primce factors...
    _NewestCheckNumbers_by_cover_Subroutine(qf_data, T_list, num_of_prime_factors, possible_exception_list_vec[XXXXX]);  
    */

  }  


  // Otherwise, use the Quadratic exceptions routine to finish them off.
  else {

    // Concatenate all possible exceptions into one list
    vector<mpz_class> possible_exception_list;
    for(long ind=1; ind < _max_num_prime_factors; ind++)
      possible_exception_list.insert(possible_exception_list.end(), possible_exception_list_vec[ind].begin(), possible_exception_list_vec[ind].end());

    
    // Status report after that checking...
    cout << endl;
    cout << endl;
    cout << "The statistics for checking with this local cover (#" << cover_number << ") are:" << endl;
    cout << endl;
    cout << _local_stats_vec[cover_number-1] << endl;    // Note: We start indexing from 0 here.
    cout << endl;
    cout << endl;
    cout << " There are " << possible_exception_list.size() << " remaining possible exceptions." << endl;
    cout << "   They are: " << possible_exception_list << endl;
    cout << endl;
    cout << " The largest possible exception " << max_eligible << " is at most " << BIGGEST_EXCEPTION_FOR_EASY_CHECKING << endl;
    cout << " so we're going to check these with the Quaternary Exceptions class! =)" << endl;
    cout << endl;
   

    // Check the (small number of) possible exceptions using the QuaternaryExceptions class! =)
    set<mpz_class> small_set;
    for(long ind=0; ind < possible_exception_list.size(); ind++) 
      small_set.insert(possible_exception_list[ind]);
    
    set<mpz_class> done_set;
    QuaternaryExceptions SMALL_EXCEPTIONS(_BigForm, small_set);
    done_set = SMALL_EXCEPTIONS.GetExceptions();

    // Sanity Check
    assert(_squarefree_exceptions.empty() == true);      

    for(set<mpz_class>::iterator ind = done_set.begin(); ind != done_set.end(); ind++)
      _squarefree_exceptions.push_back(*ind);  
    
  }










  /*
  // Write the list of exceptions...
  cout << "There are " << exception_list.size() << " exceptions. " << endl;
  cout << "They are: " << endl;
  cout << "[ ";
  if (exception_list.size() > 0) {
    for(unsigned long j=0; j < (exception_list.size()-1); j++) 
      cout << exception_list[j] << ", ";
    cout << exception_list[exception_list.size()-1];
  }
  cout << "] " << endl;


  // Write the list of overflows...
  cout << "There are " << overflow_list.size() << " overflows. " << endl;
  cout << "They are: " << endl;
  cout << "[ ";
  if (overflow_list.size() > 0) {
    for(unsigned long j=0; j < (overflow_list.size()-1); j++) 
      cout << overflow_list[j] << ", ";
    cout << overflow_list[overflow_list.size()-1];
  }
  cout << "] " << endl;
  */

  // Sort and Return the exceptions
  sort(_squarefree_exceptions.begin(), _squarefree_exceptions.end());  


  // Write the (square-free) exceptions if we're supposed to.
  _WriteSquarefree(qf_data);


  // ------ Find the square-free exceptions of the original form (not just the sublattice) --------

  // Compute the theta function of Q up to the largest exception
  //   This is the last exception:  exception_list[exception_list.size()-1]


  
  // Write the exceptions (for both) to files
  

}






// Check the representability of a vector of eligible numbers.
void Representability::_NewestCheckNumbers_by_cover_Subroutine(const QF_Datafiles & qf_data, const vector<mpz_class> & num, 
							       const long & num_of_prime_factors, vector<mpz_class> & possible_exception_list, 
							       vector<mpz_class> & possible_exception_depths) { 


  /*
  // Clear the vector to return
  possible_exception_list.clear();     // WARNING: THIS IS A BAD IDEA' since we need to keep track of these!
  */

  // /*
  // DIAGNOSTIC -- special_M
  string special_M_str("35923793310615210");
  cout << " special_M_str = " << special_M_str << endl; 
  //mpz_class special_M(290);
  mpz_class special_M(special_M_str.c_str());
  // */

  
  // Check each of the numbers in num
  for (unsigned long i=0; i< num.size(); i++) {

    bool done_flag = false;                      // This says when to finish using ternaries to check representability of a number
    bool repn_flag = false;                      // This says when a number is represented by some split local covering form
    unsigned long cover_index = _current_cover_num - 1;   
    unsigned long ternary_form_index = 0;
    //    const long MAXIMUM_NUMBER_OF_LOCAL_COVERS = 10;          // This is the maximum of local covers we'll use before giving up! =) 
    const long BIGGEST_EXCEPTION = 50000;                    // This is the bound for how large an exception must be to look for another local cover!

    long depth;               // Keeps track of the depth at which we check an eligible number

    // /*
    // DIAGNOSTIC -- special_M
    if (num[i] == special_M)
      cout << " STARTING TO CHECK m = " << special_M << endl;
    // */


    
    // STATS: Count this eligible number in the statistics for this cover
    _local_stats_vec[cover_index].total_eligible[num_of_prime_factors]++;


    // Run through all ternaries in our current cover, checking representability 
    while ((repn_flag == false) && (ternary_form_index < _current_cover_ternary_thetas.size())) {
      
      bool overflow_flag = false;  // This detects an overflow of the difference m - dx^2 for the ternary form      
      depth = 0;
      
      // Then run through all allowed depths, checking the differences for representability  
      while ((repn_flag == false) && (overflow_flag == false)) {                                  // FIX THIS TO ALLOW ANY DEPTH UNTIL OVERFLOW!!! 
	
	// Find the w to use when subtracting d*(w)^2 
	mpz_class w;
	long d = _local_cover.Get_Value(cover_index, ternary_form_index);
	w = mpz_class(floor(sqrt(num[i].get_d() / d))) - depth;    // Check this is correct...
	if (w < 0) {
	  w = 0;	
	  overflow_flag = true;    // This is really an underflow, but it still means we should stop checking!
	}

	// Find the difference
	mpz_class diff;
	diff = num[i] - mpz_class(d) * w * w;
	
	// /*
	// DIAGNOSTIC -- special_M
	if (num[i] == special_M) {
	cout << " ternary form index = " << ternary_form_index << endl;
	cout << " ternary_theta.size() = " << _ternary_theta.size() << endl;
	cout << " depth = " << depth << endl;
	cout << "   d = " << d << endl;
	cout << "   i = " << i 
	<< "   num[i] = " << num[i] << endl;
	cout << "   w = " << w 
	<< "   diff =  " << diff << endl;
	cout << endl;
	}
	// */
	
	
	// Sanity Checks
	if (diff < 0) {
	  cout << "Error in computing the difference for num[" << i <<"] = "<< num[i] << endl;
	  exit(1);
	}
	
	
	
	// Check if the difference is represented by the current ternary:
	// --------------------------------------------------------------
	
	// Check that each difference is in range
	if (diff <= _current_cover_ternary_thetas[ternary_form_index].precision()) {
	  
	  // If so, check if the ternary theta series represents the difference
	  if (_current_cover_ternary_thetas[ternary_form_index].get_value(diff) == true) {
	    repn_flag = true;

	    /*
	    // DIAGNOSTIC 
	    cout << num[i] << " is REPRESENTED!" << endl;
	    cout << " cover_index = " << cover_index << endl;
	    cout << " num_of_prime_factors = " << num_of_prime_factors << endl;
	    cout << " depth = " << depth << endl;
	    cout << " _local_stats_vec.size() = " 
		 << _local_stats_vec.size() << endl;
	    cout << " _local_stats_vec[cover_index].representation_depths.size() = " 
		 << _local_stats_vec[cover_index].representation_depths.size() << endl;
	    cout << " _local_stats_vec[cover_index].representation_depths[num_of_prime_factors].size() = " 
		 << _local_stats_vec[cover_index].representation_depths[num_of_prime_factors].size() << endl;
	    */
	    
	    
	    if (depth <= _local_stats_vec[cover_index].MAX_STATS_DEPTH)
	      _local_stats_vec[cover_index].representation_depths[num_of_prime_factors][depth]++;  // STATS: Count this number as represented at this depth   
	    
	  }
	  
	  
	  // cout << " ------ Looking at the number " << num[i] << endl;
	  
	  // /*
	  // DIAGNOSTIC -- special_M
	  if (num[i] == special_M) {
	  cout << " --> repn_flag = " << repn_flag << endl;
	  
	  cout << endl;
	  for (unsigned long mm=0; mm<=30; mm++)
	  cout << mm << " -> " << _current_cover_ternary_thetas[ternary_form_index].get_value(mm) << endl;
	  cout  << endl;
	  }
	  // */


	  // Finished checking with this depth, so increment
	  depth++;
	  
	}
	else 
	  overflow_flag = true;   // Say we overflow the precision of this form...

      }
      
      // Finished checking with this ternary form, so increment
      ternary_form_index++;
    }
    
    

    /*
      
      // Deal with being missed by all existing local covers
      if (repn_flag == false) {

	// Make another local cover if we missed a big number, but ran out of local covers.
	if ((num[i] >= BIGGEST_EXCEPTION) && (_ternary_theta.size() < MAXIMUM_NUMBER_OF_LOCAL_COVERS)) {
	  cout << " ... we missed the number " << num[i] << " with our existing local covers..." << endl;
	  
	  // Add another local cover
	  cout << " Adding another local cover... ";
	  _local_cover.IncrementLocalSplitCoverings();      
	  cout << "  giving a total of " << _local_cover.multiplicity() << " split local coverings " << endl << endl;
	  
	  // Make the ternary theta function for the maximum allowed depth
	  cout << " Computing the relevant boolean ternary theta series for cover #" << cover_number << " : " << endl;
	  cout << " ------------------------------------------------------------------- " << endl;
	  _UpdateTernaryThetas(qf_data);
	  cout << " Finished computing the relevant boolean ternary theta series for this iteration " << endl << endl << endl;
	  
	}

	// Otherwise, exit with this as a possible exception!
	else
	  done_flag = true;

      } 
      

      */
    
    
    
    
    
    
    // /*
    // DIAGNOSTIC -- special_M
    if (num[i] == special_M) {
    cout << " repn_flag = " << repn_flag << endl;
    cout << " Is it locally represented? " << _local_cover.Big_Local().IsLocallyRepresented(num[i]) << endl;
    cout << endl;
    }
    // */
    

    /*    
    // SANITY CHECK -- we should only see locally represetned numbers!
    assert(_local_cover.Big_Local().IsLocallyRepresented(num[i]) == true);
    */

    
    // If the number is not represented using any ternary (and not an overflow), 
    // but it is locally represented, then add it to the exception list. =)
    if (repn_flag == false) {
      possible_exception_list.push_back(num[i]);	     
      possible_exception_depths.push_back(depth);
      _local_stats_vec[cover_index].total_missed[num_of_prime_factors]++;        // STATS: Count this number as missed    
    }
    else 
      _local_stats_vec[cover_index].total_represented[num_of_prime_factors]++;   // STATS: Count this number as represented    
    

  }
  
  
}






/*

// Check the representability of a vector of eligible numbers by
// increasing the ternary precision until the possible exceptions are
// all small, then using the theta function to check them.
void Representability::_NewestCheckNumbers_by_size(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const int Max_Depth) {


  // Here is where we deal with overflows and decide which possible 
  // square-free exceptions and overflows are genuine exceptions!
  // --------------------------------------------------------------


  // Declare some constants and local variables
  vector<mpz_class> possible_exception_list, overflow_list; 
  const long BIGGEST_EXCEPTION = _BIGFORM_THETA_PRECISION;       // Reduce the possible exceptions to be <= 10 thousand
  vector<mpz_class> new_eligible_list;
  mpz_class max_eligible;
  long iter_depth = 1;
  bool ternary_precision_ok = true;  // Says when the ternary precision has gotten too large!

  
  // Continue checking so long as the largest candidate for an exception is too big
  while ((iter_depth == 1) || ((new_eligible_list.empty() == false) && (max_eligible > BIGGEST_EXCEPTION) && (ternary_precision_ok == true))) {   

    // Make the relevant boolean ternary theta series
    cout << " Computing the relevant boolean ternary theta series for iteration #" << iter_depth << " : " << endl;
    cout << " ----------------------------------------------------------------------- " << endl;
    ternary_precision_ok = _MakeTernaryThetas(qf_data, iter_depth);    
    cout << " Finished computing the relevant boolean ternary theta series for this iteration " << endl << endl << endl;

    // Check the representability of the numbers in num using the (iteration_depth=1) ternary theta functions
    if (iter_depth == 1)
      _NewestCheckNumbers_Subroutine(num, Max_Depth, possible_exception_list, overflow_list);
    else 
      _NewestCheckNumbers_Subroutine(new_eligible_list, Max_Depth, possible_exception_list, overflow_list);


    // List the eligible exceptions  (from the possible exceptions and overflows)
    new_eligible_list = possible_exception_list;
    new_eligible_list.insert(new_eligible_list.end(), overflow_list.begin(), overflow_list.end());
    sort(new_eligible_list.begin(), new_eligible_list.end());
    
    // Find the largest eligible exception (if there are some)
    if (new_eligible_list.empty() == false)
      max_eligible = new_eligible_list[new_eligible_list.size() - 1];


    // Increase the depth of the ternary precisions
    iter_depth++;

  }


  // Check if we've exceeded the ternary precision
  if (ternary_precision_ok == false) {

    // Update the Representability lists
    _squarefree_exceptions.insert(_squarefree_exceptions.end(), possible_exception_list.begin(), possible_exception_list.end());
    _squarefree_overflows.insert(_squarefree_overflows.end(), overflow_list.begin(), overflow_list.end());  // ??? Do we want to do this???
    
    // A not so graceful exit...
    cout << endl << endl;
    cout << " Ternary Precision overflow... Aborting..." << endl;
    cout << "   Input vector: " << num << endl;
    cout << "   Eligible Squarefree vector: " << possible_exception_list << endl;
    cout << "   Overflow vector: " << overflow_list << endl;
    cout << endl << endl;
    assert(0==1);
  }



  // Check the possible exceptions against the quaternary theta function, and append the actual exceptions
  for(long ind=0; ind < new_eligible_list.size(); ind++) 
    if (_big_theta[new_eligible_list[ind].get_ui()] == 0)
      _squarefree_exceptions.push_back(new_eligible_list[ind]);  



}

*/




// Check the representability of a vector of eligible numbers by
// increasing the ternary precision until we are reduced to a few
// small possible exceptions, then checking them individually.
void Representability::_NewestCheckNumbers_by_cover(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const long & num_of_prime_factors) {


  // Here is where we deal with overflows and decide which possible 
  // square-free exceptions and overflows are genuine exceptions!
  // --------------------------------------------------------------


  // Declare some constants and local variables
  vector<mpz_class> possible_exception_list, possible_exception_depths, overflow_list; 
  const long BIGGEST_EXCEPTION = 10000;
  const long EXCEPTION_LIMIT = 1000000;              // Don't continue if there are exceptions this big at the end.
  const long BIGGEST_NUMBER_OF_EXCEPTIONS = 100;
  vector<mpz_class> new_eligible_list;
  mpz_class max_eligible;
  long iter_depth = 1;
  bool ternary_precision_ok = true;  // Says when the ternary precision has gotten too large!

  // Check the vector of numbers
  _NewestCheckNumbers_by_cover_Subroutine(qf_data, num, num_of_prime_factors, possible_exception_list, possible_exception_depths);




  /*  
  // Continue checking so long as the largest candidate for an exception is too big
  while ((iter_depth == 1) || ((new_eligible_list.empty() == false) 
			       && (max_eligible > BIGGEST_EXCEPTION)
			       && (ternary_precision_ok == true))) {   

    // Make the relevant boolean ternary theta series
    cout << " Computing the relevant boolean ternary theta series for iteration #" << iter_depth << " : " << endl;
    cout << " ----------------------------------------------------------------------- " << endl;
    ternary_precision_ok = _MakeTernaryThetas(qf_data, iter_depth);    
    cout << " Finished computing the relevant boolean ternary theta series for this iteration " << endl << endl << endl;

    // Check the representability of the numbers in num using the (iteration_depth=1) ternary theta functions
    if (iter_depth == 1)

    else 
      _NewestCheckNumbers_Subroutine(new_eligible_list, Max_Depth, possible_exception_list);


    // List the eligible exceptions  (from the possible exceptions)
    new_eligible_list = possible_exception_list;
    sort(new_eligible_list.begin(), new_eligible_list.end());
    
    // Find the largest eligible exception (if there are some)
    if (new_eligible_list.empty() == false)
      max_eligible = new_eligible_list[new_eligible_list.size() - 1];


    // Increase the depth of the ternary precisions
    iter_depth++;

  }


  // Check if we've exceeded the ternary precision and have too many exceptions
  if ((ternary_precision_ok == false)  && ((new_eligible_list.size() > BIGGEST_NUMBER_OF_EXCEPTIONS) || (max_eligible > EXCEPTION_LIMIT)))
    {

    // Update the Representability lists
    _squarefree_exceptions.insert(_squarefree_exceptions.end(), possible_exception_list.begin(), possible_exception_list.end());
    
    // A not so graceful exit...
    cout << endl << endl;
    cout << " Ternary Precision overflow... Aborting..." << endl;
    cout << " There are " << possible_exception_list.size() << " > " << BIGGEST_NUMBER_OF_EXCEPTIONS << " remaining possible exceptions." << endl;
    //       cout << "   Input vector: " << num << endl;
    cout << "   Eligible Squarefree vector: " << possible_exception_list << endl;
    cout << "   Overflow vector: " << overflow_list << endl;
    cout << endl << endl;
    assert(0==1);
  }

  */


  // Status report after that checking...
  cout << endl;
  cout << " There are " << possible_exception_list.size() << " remaining possible exceptions." << endl;
  cout << "   They are: " << possible_exception_list << endl;
  cout << endl;





  // Check the (small number of) possible exceptions using the QuaternaryExceptions class! =)
  set<mpz_class> small_set;
  for(long ind=0; ind < possible_exception_list.size(); ind++) 
    small_set.insert(possible_exception_list[ind]);

  set<mpz_class> done_set;
  QuaternaryExceptions SMALL_EXCEPTIONS(_BigForm, small_set);
  done_set = SMALL_EXCEPTIONS.GetExceptions();

  for(set<mpz_class>::iterator ind = done_set.begin(); ind != done_set.end(); ind++)
    _squarefree_exceptions.push_back(*ind);  



}








/////////////////////////////////////////////////////////////////////////////////////////////
// Make the vector of boolean ternary theta functions for the ternaries in our local cover
/////////////////////////////////////////////////////////////////////////////////////////////

void Representability::_Get_Cover_Ternary_Thetas(const QF_Datafiles & qf_data, const long & cover_num) {
  
  // To Do:
  // ------
  // Make an ordered subset of the split local cover to use to check
  // the representability of squarefree numbers, with the best ones first! =)
  //
  // Make the option to read in the relevant boolean_ternary_thetas,
  // looking for those with greater precision, as well as making sure
  // they are in reduced form (so we don't keep any duplicates).
  //
  // These can also be the approximate theta functions... which are
  // computed much faster! =)

  /*
  // Some sanity checks
  assert(_local_cover.multiplicity() >= _ternary_theta.size());    // Check there are more local covers than computed theta functions.
  for(long i=0; i < _ternary_theta.size(); i++) {             // Check that the number of computed theta functions agrees 
    cout << " i = " << i << endl;
    cout << " _local_cover.size(i) = " << _local_cover.size(i) << endl;
    cout << " _ternary_theta[i].size() = " << _ternary_theta[i].size() << endl;
    assert(_local_cover.size(i) == _ternary_theta[i].size());      //   with the number of forms in the local covers.
  }
  */


  // Sanity Check -- Check the cover we want is in range
  assert(1 <= cover_num);  
  assert(cover_num <= _local_cover.multiplicity()); 


  //  /*
  // DIAGNOSTIC:
  cout << endl << " Entering _Get_Cover_Ternary_Thetas: " << endl;
  cout << "   Fetching the boolean theta function for cover " << cover_num << endl;
  cout << "   _ternary_theta.size() = " << _ternary_theta.size() << endl;
  cout << "   _local_cover.multiplicity() = " << _local_cover.multiplicity() << endl;
  cout << endl;
  //  */


  // Set the current local cover
  _current_cover_num = cover_num;

  // Clear the current local cover
  _current_cover_ternary_thetas.clear();

  
  
  // Set the local cover index (starting from 0)
  long i = cover_num - 1;
  
  
  // Loop through each form in this local cover
  for(long j=0; j < _local_cover.size(i); j++) {
    
    // Compute the precision for the ternaries for each local cover
    mpz_class Precision = _ThetaPrecision(_local_cover.Get_Value(i,j));      
    cout << "   using ternary precision " << Precision << endl;
    
    // Convert the current ternary form into a vector for the boolean_theta class  // <------------- *** FIX THIS RIDICULOUS NAMING SCHEME!!! ***
    long Ternary[6];
    Ternary[0] = _local_cover.Get_Form(i,j)(1,1).get_si() / 2;
    Ternary[1] = _local_cover.Get_Form(i,j)(1,2).get_si();
    Ternary[2] = _local_cover.Get_Form(i,j)(1,3).get_si();
    Ternary[3] = _local_cover.Get_Form(i,j)(2,2).get_si() / 2;
    Ternary[4] = _local_cover.Get_Form(i,j)(2,3).get_si();
    Ternary[5] = _local_cover.Get_Form(i,j)(3,3).get_si() / 2;
    
    
    
    // Declare the theta function with the necessary precision
    boolean_theta temp_theta(_local_cover.Get_Form(i,j), Precision);
    
    
    // Make the boolean theta function filename "DIR/ternary_bool_theta__a_d_f_e_c_b__precision.bindata"
    string theta_filename;
    theta_filename =  qf_data.Approx_Boolean_Theta_Dir + "ternary_bool_theta__" +
      MakeString(Ternary[0]) + "_" + MakeString(Ternary[3]) + "_" + 
      MakeString(Ternary[5]) + "_" + MakeString(Ternary[4]) + "_" + 
      MakeString(Ternary[2]) + "_" + MakeString(Ternary[1]) + 
      "__" + MakeString(Precision) + ".bindata";
    
    // DIAGNOSTIC
    cout << endl << " Looking for the theta function file: " << theta_filename << endl << endl;
    
    
    // Read the theta function, and if that fails then compute and write it!
    if ( temp_theta.read(theta_filename.c_str()) == false ) {      
      cout << endl << " Can't read the boolean theta function, so we're recomputing it!" << endl;


      // DIAGNOSTIC
      cout << " Computing the Approximate boolean theta function of " << _local_cover.Get_Form(i,j)
	   << " to precision " << Precision << endl << endl;
      
      // Use our approximate (PARI-based) routine to make the boolean theta function directly
      PrintTime();
      //temp_theta = Theta_PARI_2_new_Approximate(_local_cover.Get_Form(i,j), Precision);
      //temp_theta = Theta_PARI_3_new_Approximate_Ternary(_local_cover.Get_Form(i,j), Precision);
      temp_theta = Theta_PARI_3_new_Approximate_Ternary_mpz(_local_cover.Get_Form(i,j), Precision, APPROXIMATE_THETA_BOX_SIZE);
      PrintTime();


      // DIAGNOSTIC -- Print the first 100 numbers represented by this ternary
      vector<mpz_class> tmp_theta_numbers;
      long Test_Max = 100;
      for(long num = 0; num <= Test_Max; num++) {
	if (temp_theta.get_value(num) == 1)
	  tmp_theta_numbers.push_back(mpz_class(num));
      }
      cout << " The numbers <" << Test_Max << " are:" << endl;
      cout << tmp_theta_numbers << endl << endl;

      
      /*	
      // Get the ternary theta function from Magma
      PowerSeries<mpz_class> magma_power_series;
      magma_power_series = _local_cover.Get_Form(i,j).GetMagmaThetaSeries("", Precision);
      
      
      // Declare the theta function with the necessary precision
      boolean_ternary_theta temp_magma_theta(Ternary, magma_power_series);
      temp_theta = temp_magma_theta;
      
      // THIS IS NOT WOKING PROPERLY... =(
      //      temp_theta.compute();	
      
      */
      
      temp_theta.write(theta_filename.c_str());
    }
    
    // Add it to the temporary list of boolean ternary theta functions for this cover
    _current_cover_ternary_thetas.push_back(temp_theta);
    
  }
  
  
  
  // (Frivolous) Sanity check to ensure we have as many boolean theta functions as forms in our cover.
  assert (_current_cover_ternary_thetas.size() == _local_cover.size(i));
  
  
}

