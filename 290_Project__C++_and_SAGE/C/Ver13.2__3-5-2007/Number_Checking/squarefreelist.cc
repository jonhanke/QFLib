

//////////////////////////////////////////////////////////////////////
// Original routine to make the list of eligible squarefree numbers //
//////////////////////////////////////////////////////////////////////

void SharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP) {

  // Find out the maximal number of prime factors
  int num = 0;
  double temp_F = 1;
  while (temp_F <= bound) {
    temp_F = temp_F * F4PP[num];
    num = num + 1;
    cout << " The product after " << num << " factors is: " << temp_F << endl;
  }

  int Max_num = num - 1;
  cout <<  "We can have square-free numbers with at most " << Max_num << " prime factors." << endl;


  // Make the array of Max_num pointers and step through them 
  // to get a minimal list of square-free numbers t to check.

  long PP_length = PP.size();
  int carry_ptr;
  int counter=0;
  

  // Loop through all the possible #'s of prime factors
  for (int num = Max_num; num >= 1; num--) {
    
    // Create the output filename: 
    //    "eligible_square-free__bound_number-of-primes__.data"
    char T_filename[70];
    sprintf(T_filename, "eligible_square-free__b=%f_total-p=%lu_factors=%d___.data", bound, PP.size(),num); 
    //cout << "\n Filename: " << T_filename << "X" << endl << endl;  
        
    // Open the output file and check it's opened correctly
    ofstream T_file; 
    T_file.open(T_filename);
    
    if (! T_file.is_open())
      { cout << "Error opening output file " << T_filename; exit (1); }

    
    // Initialize the array for j prime factors
    int Last_ptr = num - 1;
    carry_ptr = Last_ptr;  // This keeps the compiler from complaining... =)
    vector<long> Ptr_array(num);
    for (int i=0; i<= Last_ptr; i++) 
      Ptr_array[i] = i;

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
	  T_file << temp_squarefree << ", ";
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
    }


    // /*
    // Description of the number of square-free numbers with num prime factors:
    cout << "\n Decreasing the number of prime factors." << endl;
    cout << "   Last_ptr being decreased from " << Last_ptr << " to " << Last_ptr - 1 << endl;
    cout << "   Added " << T_list.size() - counter << " square-free numbers with " 
	 << num << " prime factors." << endl;
    cout << "   The T_list has " << T_list.size() << " numbers." << endl;
    PrintTime();

    counter = T_list.size();  // Reset the counter for the next round. 
    // */


    // Remove the trailing ", " and insert an EOF
    long file_ptr;
    file_ptr = T_file.tellp();
    T_file.seekp(file_ptr - 2);
    T_file.put(' ');
    T_file.put(EOF);
    T_file.close();

  }
  
  //  return T_list;
}





////////////////////////////////////////////////////////////////
// Forgetful version which doesn't keep the list in memory =) //
////////////////////////////////////////////////////////////////

void ForgetfulSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP) {

  // Find out the maximal number of prime factors
  int num = 0;
  double temp_F = 1;
  while (temp_F <= bound) {
    temp_F = temp_F * F4PP[num];
    num = num + 1;
    cout << " The product after " << num << " factors is: " << temp_F << endl;
  }

  int Max_num = num - 1;
  cout <<  "We can have square-free numbers with at most " << Max_num << " prime factors." << endl;


  // Make the array of Max_num pointers and step through them 
  // to get a minimal list of square-free numbers t to check.

  long PP_length = PP.size();
  int carry_ptr;
  int counter=0;
  

  // Loop through all the possible #'s of prime factors
  for (int num = Max_num; num >= 1; num--) {
    
    // Create the output filename: 
    //    "eligible_square-free__bound_number-of-primes__.data"
    char T_filename[70];
    sprintf(T_filename, "eligible_square-free__b=%f_total-p=%lu_factors=%d___.data", bound, PP.size(),num); 
    //cout << "\n Filename: " << T_filename << "X" << endl << endl;  
        
    // Open the output file and check it's opened correctly
    ofstream T_file; 
    T_file.open(T_filename);
    
    if (! T_file.is_open())
      { cout << "Error opening output file " << T_filename; exit (1); }

    
    // Initialize the array for j prime factors
    int Last_ptr = num - 1;
    carry_ptr = Last_ptr;  // This keeps the compiler from complaining... =)
    vector<long> Ptr_array(num);
    for (int i=0; i<= Last_ptr; i++) 
      Ptr_array[i] = i;

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
                    
	  //T_list.push_back(temp_squarefree);
	  T_file << temp_squarefree << ", ";
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
    }


    // /*
    // Description of the number of square-free numbers with num prime factors:
    cout << "\n Decreasing the number of prime factors." << endl;
    cout << "   Last_ptr being decreased from " << Last_ptr << " to " << Last_ptr - 1 << endl;
    cout << "   Added " << T_list.size() - counter << " square-free numbers with " 
	 << num << " prime factors." << endl;
    cout << "   The T_list has " << T_list.size() << " numbers." << endl;
    PrintTime();

    counter = T_list.size();  // Reset the counter for the next round. 
    // */


    // Remove the trailing ", " and insert an EOF
    long file_ptr;
    file_ptr = T_file.tellp();
    T_file.seekp(file_ptr - 2);
    T_file.put(' ');
    T_file.put(EOF);
    T_file.close();

  }
  
  //  return T_list;
}



///////////////////////////////////////////////////////////////
// Lazy version which doesn't write the list to the disk. =) //
///////////////////////////////////////////////////////////////

void LazySharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP) {

  // Find out the maximal number of prime factors
  int num = 0;
  double temp_F = 1;
  while (temp_F <= bound) {
    temp_F = temp_F * F4PP[num];
    num = num + 1;
    cout << " The product after " << num << " factors is: " << temp_F << endl;
  }

  int Max_num = num - 1;
  cout <<  "We can have square-free numbers with at most " << Max_num << " prime factors." << endl;


  // Make the array of Max_num pointers and step through them 
  // to get a minimal list of square-free numbers t to check.

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
    }


    // /*
    // Description of the number of square-free numbers with num prime factors:
    cout << "\n Decreasing the number of prime factors." << endl;
    cout << "   Last_ptr being decreased from " << Last_ptr << " to " << Last_ptr - 1 << endl;
    cout << "   Added " << T_list.size() - counter << " square-free numbers with " 
	 << num << " prime factors." << endl;
    cout << "   The T_list has " << T_list.size() << " numbers." << endl;
    PrintTime();

    counter = T_list.size();  // Reset the counter for the next round. 
    // */


  }
  
  //  return T_list;
}






////////////////////////////////////////////////////////////////
// An incremental version which checks as it goes along... =) //
////////////////////////////////////////////////////////////////

void NewSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP, 
		  vector<bool> & Ternary_series, long diag_coeff) {

  // Make a counter to count squarefree numbers
  // and lists of exceptions and overflows...
  mpz_class Counter;
  Counter = 0;
  vector<mpz_class> exception_list, overflow_list, overflow_differences, overflow_depths;


  // Find out the maximal number of prime factors
  int num = 0;
  double temp_F = 1;
  while (temp_F <= bound) {
    temp_F = temp_F * F4PP[num];
    num = num + 1;
    cout << " The product after " << num << " factors is: " << temp_F << endl;
  }

  int Max_num = num - 1;
  cout <<  "We can have square-free numbers with at most " << Max_num << " prime factors." << endl;


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
	NewCheckNumbers(T_list, Checking_Depth, Ternary_series, diag_coeff,
			exception_list, overflow_list, overflow_differences, overflow_depths);	
	T_list.resize(0);

	cout << "      Checked " << Counter << " numbers" << endl;
	cout << "      Number of exceptions: " << exception_list.size() << endl;
	cout << "      Number of overflows: " << overflow_list.size() << endl;
	cout << "      "; PrintTime();
	cout << endl;
      }

    }


    // /*
    // Description of the number of square-free numbers with num prime factors:
    cout << "\n Decreasing the number of prime factors." << endl;
    cout << "   Last_ptr being decreased from " << Last_ptr << " to " << Last_ptr - 1 << endl;
    cout << "   Added " << T_list.size() - counter << " square-free numbers with " 
	 << num << " prime factors." << endl;
    cout << "   The T_list has " << T_list.size() << " numbers." << endl;
    PrintTime();

    counter = T_list.size();  // Reset the counter for the next round. 
    // */


  }
  
  
  // Write the list of exceptions...
  cout << "There are " << exception_list.size() << " exceptions. " << endl;
  cout << "They are: " << endl;
  cout << "[ ";
  for(unsigned long j=0; j < (exception_list.size()-1); j++) 
    cout << exception_list[j] << ", ";
  cout << exception_list[exception_list.size()-1] << "] ";


  //  return T_list;
}






////////////////////////////////////////////////////////////////
// An incremental version which checks as it goes along... =) //
////////////////////////////////////////////////////////////////

void NewerSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP, 
		  unsigned long Ternary_series[], unsigned long Precision, long diag_coeff) {

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
	NewerCheckNumbers(T_list, Checking_Depth, Ternary_series, Precision, diag_coeff,
			   exception_list, overflow_list, overflow_differences, overflow_depths);	

	T_list.resize(0);

	cout << "      Checked " << Counter << " numbers" << endl;
	cout << "      Number of exceptions: " << exception_list.size() << endl;
	cout << "      Number of overflows: " << overflow_list.size() << endl;
	cout << "      "; PrintTime();
	cout << endl;
      }
	

    }


    // /*
    // Description of the number of square-free numbers with num prime factors:
    cout << "\n Decreasing the number of prime factors." << endl;
    cout << "   Last_ptr being decreased from " << Last_ptr << " to " << Last_ptr - 1 << endl;
    cout << "   Added " << T_list.size() - counter << " square-free numbers with " 
	 << num << " prime factors." << endl;
    cout << "   The T_list has " << T_list.size() << " numbers." << endl;
    PrintTime();

    counter = T_list.size();  // Reset the counter for the next round. 
    // */


  }
  
  
  // Write the list of exceptions...
  cout << "There are " << exception_list.size() << " exceptions. " << endl;
  cout << "They are: " << endl;
  cout << "[ ";
  if (exception_list.size() > 0) {
    for(unsigned long j=0; j < (exception_list.size()-1); j++) 
      cout << exception_list[j] << ", ";
    cout << exception_list[exception_list.size()-1];
  }
  cout << "] ";


  //  return T_list;
}













// =========================================================================================================














///////////////////////////////////////////////////////////////////////////////////////
// Modified routine to check the bound overflow after multiplying each new factor... //
///////////////////////////////////////////////////////////////////////////////////////

void FasterSharpList(vector<mpz_class> & T_list, double bound, const vector<long> & PP, const vector<double> & F4PP) {

  // Find out the maximal number of prime factors
  int num = 0;
  double temp_F = 1;
  while (temp_F <= bound) {
    temp_F = temp_F * F4PP[num];
    num = num + 1;
    cout << " The product after " << num << " factors is: " << temp_F << endl;
  }

  int Max_num = num - 1;
  cout <<  "We can have square-free numbers with at most " << Max_num << " prime factors." << endl;


  // Make the array of Max_num pointers and step through them 
  // to get a minimal list of square-free numbers t to check.

  long PP_length = PP.size();
  int carry_ptr;
  int counter=0;
  

  // Loop through all the possible #'s of prime factors
  for (int num = Max_num; num >= 1; num--) {
    
    // Create the output filename: 
    //    "eligible_square-free__bound_number-of-primes__.data1"
    char T_filename[70];
    sprintf(T_filename, "eligible_square-free__b=%f_total-p=%lu_factors=%d___.data1", bound, PP.size(),num); 
    //cout << "\n Filename: " << T_filename << "X" << endl << endl;  
        
    // Open the output file and check it's opened correctly
    ofstream T_file; 
    T_file.open(T_filename);
    
    if (! T_file.is_open())
      { cout << "Error opening output file " << T_filename; exit (1); }

    
    // Initialize the array for j prime factors
    int Last_ptr = num - 1;
    carry_ptr = Last_ptr;  // This keeps the compiler from complaining... =)
    vector<long> Ptr_array(num);
    for (int i=0; i<= Last_ptr; i++) 
      Ptr_array[i] = i;

    bool done_flag = false;
 
    while (done_flag == false) {

      bool carry_flag = false;    // Start without wanting to carry.

      // Check that we're still in the allowed primes range
      if (Ptr_array[Last_ptr] <= PP_length - 1) { 
            
	// Compute F4(t) for the product of the currently indexed primes
	// (unless there is an overflow...)
	double temp_F4_prod = 1;
	int jj = 0;
	while ((temp_F4_prod <= bound) && (jj <= Last_ptr)) {
	  temp_F4_prod = temp_F4_prod * F4PP[Ptr_array[jj]];  
	  jj = jj + 1;	  
	}

	// If F4(t) < bound, add t to the list and reset the carry_ptr 
	if (temp_F4_prod <= bound) {
	  mpz_class temp_squarefree = 1;
	  for (int i=0; i<=Last_ptr; i++) 
	    temp_squarefree = temp_squarefree * PP[Ptr_array[i]];
                    
	  T_list.push_back(temp_squarefree);
	  T_file << temp_squarefree << ", ";
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
    }


    // /*
    // Description of the number of square-free numbers with num prime factors:
    cout << "\n Decreasing the number of prime factors." << endl;
    cout << "   Last_ptr being decreased from " << Last_ptr << " to " << Last_ptr - 1 << endl;
    cout << "   Added " << T_list.size() - counter << " square-free numbers with " 
	 << num << " prime factors." << endl;
    cout << "   The T_list has " << T_list.size() << " numbers." << endl;
    PrintTime();

    counter = T_list.size();  // Reset the counter for the next round. 
    // */


    // Remove the trailing ", " and insert an EOF
    long file_ptr;
    file_ptr = T_file.tellp();
    T_file.seekp(file_ptr - 2);
    T_file.put(' ');
    T_file.put(EOF);
    T_file.close();

  }
  
  //  return T_list;
}

