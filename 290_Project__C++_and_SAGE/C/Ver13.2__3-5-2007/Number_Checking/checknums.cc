

// Check the representability of a vector of eligible numbers.
vector<mpz_class> CheckNumbers(const vector<mpz_class> & num, int depth, const vector<bool> & theta_vec, long diag_coeff) {
  
  double ddiag_coeff = double(diag_coeff);
  mpz_class zdiag_coeff = mpz_class(diag_coeff);

  // Make the list of eligible indices to check the representability of 
  // (start with all if them!)
  vector<unsigned long>  num_index_list(num.size());
  for (unsigned long i=0; i< num.size(); i++) 
    num_index_list[i] = i;

  
  // Loop through differences of depth d, checking the eligible indices for representability
  for(long d=0; d < depth; d++) { 

    // Report the status
    cout << " Checking the representability of numbers with depth " << d << endl;
    cout << " There are " << num_index_list.size() << " numbers to check." << endl << endl;
    PrintTime();
    
    // Make the list of differences with depth d
    vector<mpz_class> Dlist(num_index_list.size()); 
    mpz_class w;
    for (unsigned long i=0; i< num_index_list.size(); i++) {
      w = mpz_class(floor(sqrt(num[num_index_list[i]].get_d() / ddiag_coeff))) - d;    // Check this is correct...
      if (w < 0)  
	w=0;	
      Dlist[i] = num[num_index_list[i]] - zdiag_coeff * w * w;
      /*
      cout << "  i = " << i 
	   << "  num_index_list[i] = " << num_index_list[i] 
	   << "  num[num_index_list[i]] = " << num[num_index_list[i]] 
	   << "  w = " << w 
	   << "  Dlist[i] = " << Dlist[i] << endl;
      */
      if (Dlist[i] < 0)
	cout << "Error in computing Dlist[" << i <<"] = "<< Dlist[i] << endl;
    }
    
    
    // Compute the maximum of the Dlist, which must be less than the theta series length...
    mpz_class Dlist_max=0;
    for (unsigned long i=0; i< Dlist.size(); i++) 
      if (Dlist[i] > Dlist_max)
	Dlist_max = Dlist[i];

    if ((Dlist_max.fits_uint_p() ) && (Dlist_max < theta_vec.size())) 
      cout << " The maximum of the Dlist is " << Dlist_max 
	   << " which is less than the theta series length "<< theta_vec.size() 
	   << " and the sizeof(long) " << endl;
    else {
      cout << " ERROR: The maximum of the Dlist is " << Dlist_max 
	   << " which is >= either the theta series length "<< theta_vec.size() 
	   << " or the sizeof(long) " << "!!!" << endl;
      exit(1);   // Abort the program!
    }

    
    // Check these against the theta series to produce a refined list of indices to check
    //    vector<long> missed_list_indices;
    vector<unsigned long>  new_indices;
    for (unsigned long i=0; i< Dlist.size(); i++) 
      if (theta_vec[Dlist[i].get_ui()] == false)    // Keep it if the theta series doesn't represent the difference
	new_indices.push_back(num_index_list[i]);        
    num_index_list = new_indices;

  }

  // Make the list of numbers that weren't found to be represented 
  vector<mpz_class> failed_numbers(num_index_list.size());
  for (unsigned long i=0; i< num_index_list.size(); i++) 
    failed_numbers[i] = num[num_index_list[i]];
  
  return failed_numbers;  
}






// Check the representability of a vector of eligible numbers.
void NewCheckNumbers(const vector<mpz_class> & num, int depth, const vector<bool> & theta_vec, long diag_coeff, 
		vector<mpz_class> & exception_list, vector<mpz_class> & overflow_list, 
		vector<mpz_class> & overflow_values, vector<mpz_class> & overflow_depths) {
  
  double ddiag_coeff = double(diag_coeff);
  mpz_class zdiag_coeff = mpz_class(diag_coeff);

  // Make the list of eligible indices to check the representability of 
  // (start with all if them!)
  vector<unsigned long>  num_index_list(num.size());
  for (unsigned long i=0; i< num.size(); i++) 
    num_index_list[i] = i;

  
  // Loop through differences of depth d, checking the eligible indices for representability
  for(long d=0; d < depth; d++) { 

    // Report the status
    /*
    cout << " Checking the representability of numbers with depth " << d << endl;
    cout << " There are " << num_index_list.size() << " numbers to check." << endl << endl;
    PrintTime();
    */    

    // Make the list of differences with depth d
    vector<mpz_class> Dlist(num_index_list.size()); 
    mpz_class w;
    for (unsigned long i=0; i< num_index_list.size(); i++) {
      w = mpz_class(floor(sqrt(num[num_index_list[i]].get_d() / ddiag_coeff))) - d;    // Check this is correct...
      if (w < 0)  
	w=0;	
      Dlist[i] = num[num_index_list[i]] - zdiag_coeff * w * w;
      /*
      cout << "  i = " << i 
	   << "  num_index_list[i] = " << num_index_list[i] 
	   << "  num[num_index_list[i]] = " << num[num_index_list[i]] 
	   << "  w = " << w 
	   << "  Dlist[i] = " << Dlist[i] << endl;
      */
      if (Dlist[i] < 0)
	cout << "Error in computing Dlist[" << i <<"] = "<< Dlist[i] << endl;
    }
    
    
    // Check these against the theta series to produce a refined list of indices to check
    //    vector<long> missed_list_indices;
    vector<unsigned long>  new_indices;
    for (unsigned long i=0; i< Dlist.size(); i++) {
      
      // Check that each difference is in range
      if ((Dlist[i].fits_uint_p() ) && (Dlist[i] < theta_vec.size())) {
	
	// If so, check if the ternary theta series represents the difference
	if (theta_vec[Dlist[i].get_ui()] == false)    
	  new_indices.push_back(num_index_list[i]);        
      }

      // If not, add it to the overflow lists
      else {
	overflow_list.push_back(num[num_index_list[i]]);
	overflow_values.push_back(Dlist[i]);
	overflow_depths.push_back(d);
      }
      
    }
    num_index_list = new_indices;
    
  }

  // Append the numbers that weren't found to exception_list 
  for (unsigned long i=0; i< num_index_list.size(); i++) 
    exception_list.push_back(num[num_index_list[i]]);
  

  // Describe the results of this routine... (maybe later)

}





// Check the representability of a vector of eligible numbers.
void NewerCheckNumbers(const vector<mpz_class> & num, int depth, const unsigned long theta_vec[], 
		       unsigned long Theta_Precision, long diag_coeff, 
		       vector<mpz_class> & exception_list, vector<mpz_class> & overflow_list, 
		       vector<mpz_class> & overflow_values, vector<mpz_class> & overflow_depths) {
  
  double ddiag_coeff = double(diag_coeff);
  mpz_class zdiag_coeff = mpz_class(diag_coeff);

  

  // Make the list of eligible indices to check the representability of 
  // (start with all if them!)
  vector<unsigned long>  num_index_list(num.size());
  for (unsigned long i=0; i< num.size(); i++) 
    num_index_list[i] = i;

  
  // Loop through differences of depth d, checking the eligible indices for representability
  for(long d=0; d < depth; d++) { 

    // Report the status
    /*
    cout << " Checking the representability of numbers with depth " << d << endl;
    cout << " There are " << num_index_list.size() << " numbers to check." << endl << endl;
    PrintTime();
    */    

    // Make the list of differences with depth d
    vector<mpz_class> Dlist(num_index_list.size()); 
    mpz_class w;
    for (unsigned long i=0; i< num_index_list.size(); i++) {
      w = mpz_class(floor(sqrt(num[num_index_list[i]].get_d() / ddiag_coeff))) - d;    // Check this is correct...
      if (w < 0)  
	w=0;	
      Dlist[i] = num[num_index_list[i]] - zdiag_coeff * w * w;
      /*
      cout << "  i = " << i 
	   << "  num_index_list[i] = " << num_index_list[i] 
	   << "  num[num_index_list[i]] = " << num[num_index_list[i]] 
	   << "  w = " << w 
	   << "  Dlist[i] = " << Dlist[i] << endl;
      */
      if (Dlist[i] < 0)
	cout << "Error in computing Dlist[" << i <<"] = "<< Dlist[i] << endl;
    }
    
    
    // Check these against the theta series to produce a refined list of indices to check
    //    vector<long> missed_list_indices;
    vector<unsigned long>  new_indices;
    bool bin_bit;
    unsigned long j;
    for (unsigned long i=0; i< Dlist.size(); i++) {
      
      // Check that each difference is in range
      if ((Dlist[i].fits_uint_p() ) && (Dlist[i] <= Theta_Precision)) {
	
	// If so, check if the ternary theta series represents the difference
	j = Dlist[i].get_ui();
	bin_bit = ((theta_vec[j >> 5] & (unsigned long) (1 << (j % 32)))) >> (j % 32);
	if (bin_bit == false)    
	  new_indices.push_back(num_index_list[i]);        
      }

      // If not, add it to the overflow lists
      else {
	overflow_list.push_back(num[num_index_list[i]]);
	overflow_values.push_back(Dlist[i]);
	overflow_depths.push_back(d);
      }
      
    }
    num_index_list = new_indices;
    
  }

  // Append the numbers that weren't found to exception_list 
  for (unsigned long i=0; i< num_index_list.size(); i++) 
    exception_list.push_back(num[num_index_list[i]]);
  

  // Describe the results of this routine... (maybe later)

}





