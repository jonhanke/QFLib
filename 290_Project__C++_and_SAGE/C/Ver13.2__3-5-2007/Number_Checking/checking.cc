



//! Check the representability of a vector of eligible numbers.
void NewestCheckNumbers(const vector<mpz_class> & num, int depth, 
			boolean_ternary_theta  & theta3,
			long diag_coeff, 
			vector<mpz_class> & exception_list, vector<mpz_class> & overflow_list, 
			vector<mpz_class> & overflow_values, vector<mpz_class> & overflow_depths,
			vector<long> & local_mod_vector, long local_repn_array[200][9]) {

  
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
    unsigned long j;
    for (unsigned long i=0; i< Dlist.size(); i++) {
      
      // Check that each difference is in range
      if ((Dlist[i].fits_uint_p() ) && (Dlist[i] <= theta3.precision())) {
	
	// If so, check if the ternary theta series represents the difference
	j = Dlist[i].get_ui();
	if (theta3.get_value(j) == false)    
	  new_indices.push_back(num_index_list[i]);        
      }

      // If not, add it to the overflow lists
      else {
	cout << "ERROR: Overflow occurred at m = " << num[num_index_list[i]] 
	     << endl;
	overflow_list.push_back(num[num_index_list[i]]);
	overflow_values.push_back(Dlist[i]);
	overflow_depths.push_back(d);
      }
      
    }
    num_index_list = new_indices;
    
  }


  /*
  // DIAGNOSTIC -- all missed numbers
  cout << " Here are the exceptions before checking local representability " << endl;
  for (unsigned long i=0; i< num_index_list.size(); i++) 
    cout << num[num_index_list[i]] << ", ";
  cout << endl << endl;
  */
  
  
  // ------------------- Check if the missed numbers are locally represented ----------------------------
  vector<unsigned long> excep_index_list;
  
  for (unsigned long i=0; i< num_index_list.size(); i++) {

    // Set m to be the number to check the local representability of
    mpz_class m;
    m = num[num_index_list[i]];

   
    // Check for local representability at each eligible prime
    bool is_locally_missed = false;
    for (unsigned long k=0; k < local_mod_vector.size(); k++) {  // NOTE: Could speed this up by exiting as soon as is_missed is true
      
      long p = local_repn_array[k][0];
      //      unsigned long v = Valuation(m, p);
      
      //	  cout << " Hi -- Using m = " << m << " and p = " << p << endl;
      //	  cout << "  valuation = " << v << endl;
      
      
      // Find the associated number kk (so m = x * p^(2*kk))
      //      unsigned long kk = v / 2;  
      unsigned long kk = 0;  // kk is always 0;
      
      //      long m1 = m / LongPow(p, 2*kk);
      long m1 = m.get_ui(); // This is the squarefree number -- which is what m is by construction.
      
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
      excep_index_list.push_back(num_index_list[i]);
    
  }
  

  /*
  // DIAGNOSTIC -- missed numbers which are not locally represented 
  cout << " Here are the exceptions after checking local representability " << endl;
  for (unsigned long i=0; i< excep_index_list.size(); i++) 
    cout << num[excep_index_list[i]] << ", ";
  cout << endl << endl;
  */


  /*
  // DIAGNOSTIC -- overflow list
  cout << " Here are the overflows " << endl;
  for (unsigned long i=0; i< overflow_list.size(); i++) 
    cout << overflow_list[i] << ", ";
  cout << endl << endl;
  */

  
  // Append the numbers that weren't found to exception_list 
  for (unsigned long i=0; i< excep_index_list.size(); i++) 
    exception_list.push_back(num[excep_index_list[i]]);
  
  
  // Describe the results of this routine... (maybe later)
}












// An incremental version which checks as it goes along... =) 

//! An incremental version which checks as it goes along... =) 

/*! 
  This is the detailed comment... but I probably need to say something more!
  \f$\sqrt{ (x_2-x_1)^2 + (y_2 - y_1)^2 }\f$
  Let's try using the formula \f$\int_a^b f(x) dx \f$. =)
*/



vector<mpz_class> NewestSharpList(vector<mpz_class> & T_list, double bound, //< Bound for F_4(p).
				  const vector<long> & PP, const vector<double> & F4PP, 
				  boolean_ternary_theta & theta3,
				  long diag_coeff, Matrix_mpz QQ) {


  // Make a counter to count squarefree numbers
  // and lists of exceptions and overflows...
  mpz_class Counter;
  Counter = 0;
  vector<mpz_class> exception_list, overflow_list, overflow_differences, overflow_depths;



  // Find out the maximal number of prime factors
  int num = 0;
  double temp_F = 1;
  while ((temp_F <= bound) && (num < (int) PP.size())) {  // Note: This is needed when we only have a few small primes!
    temp_F = temp_F * F4PP[num];
    num = num + 1;
    cout << " The product after " << num << " factors is: " << temp_F << endl;
  }

  int Max_num = num - 1;
  cout <<  "We can have square-free numbers with at most " << Max_num << " prime factors." << endl;


  // Find the local conditions
  vector<long> local_mod_vector;
  long local_repn_array[200][9];  // NOTE: We should really not have a 200 here.  
                                  // This should be dynamic, and depend on the number of p|N.
  vector<long> aniso_vector;

  // Get the ternary form from theta3
  long *Ternary;
  //  long Ternary[6];
  Ternary = theta3.get_qf();
  
  Matrix_mpz QQ_sub;  // Make the quadratic form of the (4-dim'l) decomposable sublattice!
  QQ_sub.SetDims(4,4);
  QQ_sub(1,1) = 2 * diag_coeff;
  QQ_sub(2,2) = 2 * Ternary[0];
  QQ_sub(2,3) = Ternary[1];
  QQ_sub(3,2) = Ternary[1];
  QQ_sub(2,4) = Ternary[2];
  QQ_sub(4,2) = Ternary[2];
  QQ_sub(3,3) = 2 * Ternary[3];
  QQ_sub(3,4) = Ternary[4];
  QQ_sub(4,3) = Ternary[4];
  QQ_sub(4,4) = 2 * Ternary[5];
  
  cout << " QQ is:" << endl << theta3.get_qf() << endl;
  cout << " QQ_sub is:" << endl << QQ_sub << endl;

  cout << " Using Precision = " << theta3.precision() << endl;
    
  cout << "Finding the local conditions" << endl; 
  FindLocalConditions(QQ_sub, local_mod_vector, local_repn_array, aniso_vector); 
  cout << "Found the local conditions" << endl << endl << endl; 


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



  // IS THIS NECESSARY ANYMORE???
  // Modify the local_mod_vector at p=2 to exclude multiples of mod/4 or mod for each prime power.
  local_mod_vector[0] = local_mod_vector[0] / 4; 
  


  //////////////////////////////////////////////////////////////
  // Make the array of Max_num pointers and step through them //
  // to get a minimal list of square-free numbers t to check. //
  //////////////////////////////////////////////////////////////

  long PP_length = PP.size();
  int carry_ptr;
  int counter=0;

  // Loop through all the possible #'s of prime factors
  for (int num = Max_num; num >= 1; num--) {
    
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
      if (Ptr_array[Last_ptr] <= PP_length - 1) { 
            
	// compute F4(t) for the product of the currently indexed primes
	double temp_F4_prod = 1;
	for (int i=0; i<=Last_ptr; i++)
	  temp_F4_prod = temp_F4_prod * F4PP[Ptr_array[i]];  

	// If F4(t) < bound, add t to the list and reset the carry_ptr 
	if (temp_F4_prod <= bound) {
	  mpz_class temp_squarefree = 1;
	  for (int i=0; i<=Last_ptr; i++) 
	    temp_squarefree = temp_squarefree * PP[Ptr_array[i]];
                    
	  T_list.push_back(temp_squarefree);
	  //T_file << temp_squarefree << ", ";
	  Ptr_array[Last_ptr] = Ptr_array[Last_ptr] + 1;
	  carry_ptr = Last_ptr;
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
	int Checking_Depth = 5;

	/*
	// OVERFLOW DIAGNOSTIC...
	T_list.resize(0);
	T_list.push_back(mpz_class(940869930));
	*/

	NewestCheckNumbers(T_list, Checking_Depth, theta3, diag_coeff,
			   exception_list, overflow_list, overflow_differences, overflow_depths,
			   local_mod_vector, local_repn_array);	

	T_list.resize(0);

	cout << "     Cumulative results: " << endl;
	cout << "      Checked " << Counter << " numbers " << endl;
	cout << "      Number of exceptions: " << exception_list.size() << endl;
	cout << "      Number of overflows: " << overflow_list.size() << endl;
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
  sort(exception_list.begin(), exception_list.end());  
  return exception_list;



  // ------ Find the square-free exceptions of the original form (not just the sublattice) --------

  // Compute the theta function of Q up to the largest exception
  //   This is the last exception:  exception_list[exception_list.size()-1]


  
  // Write the exceptions (for both) to files
  


}
