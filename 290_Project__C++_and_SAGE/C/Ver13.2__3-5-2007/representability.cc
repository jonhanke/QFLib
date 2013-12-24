
///////////////////////////
// Constructor front-end
///////////////////////////
Representability::Representability(const Matrix_mpz & QQ, const double & cusp_const, const string & type) {

  // Sanity Check
  assert( (type == "uniform") || (type == "minimal"));

  // Important pre-initialization to detect errors
  _max_num_prime_factors = -1;

  // Initialize the Default Checking Depth and Approximate theta function box size
  DEFAULT_CHECKING_DEPTH = 9;                    // Default depth is 5
  APPROXIMATE_THETA_BOX_SIZE = 800.;             // Default box dimensions are 500.


  // DIAGNSOTIC
  cout << " Using the form QQ = " << endl << QQ << endl;
  

  QF_Datafiles qf_data("QF__Representability__temp_data");  // Make temporary data folders

  // Choose the method of checking 
  if (type == "minimal")
    (*this)._MakeWithMinimalCover(qf_data, QQ, cusp_const);
  else if (type == "uniform")
    (*this)._MakeWithUniformCover(qf_data, QQ, cusp_const);
  
  qf_data.Delete_All_Folders();                             // Delete the temporary folders

}



///////////////////////////
// Constructor front-end
///////////////////////////
Representability::Representability(const QF_Datafiles & qf_data, const Matrix_mpz & QQ, const double & cusp_const, const string & type) {

  // Sanity Check
  assert( (type == "uniform") || (type == "minimal"));

  // Important pre-initialization to detect errors
  _max_num_prime_factors = -1;

  // Initialize the Default Checking Depth and Approximate theta function box size
  DEFAULT_CHECKING_DEPTH = 9;                    // Default depth is 5
  APPROXIMATE_THETA_BOX_SIZE = 800.;             // Default box dimensions are 500.


  // DIAGNSOTIC
  cout << " Using the form QQ = " << endl << QQ << endl;


  // Choose the method of checking 
  if (type == "minimal")
    (*this)._MakeWithMinimalCover(qf_data, QQ, cusp_const);
  else if (type == "uniform")
    (*this)._MakeWithUniformCover(qf_data, QQ, cusp_const);
  

}


////////////////////////
// Secret Constructor
////////////////////////
void Representability::_MakeWithMinimalCover(const QF_Datafiles & qf_data, const Matrix_mpz & QQ, const double & cusp_const) {

  // TO DO:
  //   Check that Q is 4x4, symmetric, and positive definite


  // Set the exception size bound
  _BIGFORM_THETA_PRECISION = 10000;


  // Set the Overflow Bound
  const double OVERFLOW_BOUND = 20000;                           // ********* This should be 2000... **********
  overflow = false;    // Default value needed for when the cusp_const is zero.


  // Store the (pre-computed) cuspidal upper bound, allowing for roundoff errors! =)
  _passed_cusp_const = cusp_const + 1;
  cout << " Passed in the cusp constant " << cusp_const << endl;
  cout << " Using the larger cusp constant " << _passed_cusp_const << " to allow for roundoff errors! =)" << endl << endl;

  
  // Compute the level and the character of QQ
  _BigForm = QQ;
  _level = QQ.QFLevel().get_si();
  _character = SquarefreePart(QQ.Determinant()).get_si();
  
  // TO LOOK AT:
  //   unsigned long Precision = (unsigned long) ceil(Theta_Precisions[k]);


  // Try to read the (square-free) exceptions from a file first  <--- **** NEED TO CHANGE THIS WEHN WE CHECK SQUARE FACTORS TOO ***
  if (_ReadSquarefree(qf_data) == true) {
    _used_squarefree_exception_file = true;
    _CheckSquareFactors();
  }
  else {
    _used_squarefree_exception_file = false;
    
    
    // Skip all of the work if there are no cusp forms ! =)
    if (cusp_const != 0) {
    
      // Find the F4 Bound for QQ
      cout << " Finding the F4 Bound: " << endl;
      cout << " --------------------- " << endl;
      _GetF4Bound(qf_data);
      cout << " Finished finding the F4 Bound " << endl << endl << endl;
      
    // Check for an overflow... (when F4_bound > OVERFLOW_BOUND)
      if (_F4_bound > OVERFLOW_BOUND) {
	overflow = true;
	cout << " Representability Overflow: The bound " << _F4_bound << " is > " << OVERFLOW_BOUND << endl;
      }
      else {
	overflow = false;
	
	// Construct a local cover for QQ
	cout << " Making the minimal local cover: " << endl;
	cout << " ------------------------------- " << endl;
	_local_cover.Make("minimal", QQ, qf_data.Local_Cover_Dir);
	cout << " Finished making the minimal local cover " << endl << endl << endl;
	
	
	// Find the eligible primes 
	cout << " Finding the eligible primes: " << endl;
	cout << " ---------------------------- " << endl;
	_GetEligiblePrimes(qf_data);
	cout << " Finished finding the eligible primes " << endl << endl << endl;


	// Initialize the local cover statistics data
	cout << " Initializing the local cover statistics data: " << endl;
	cout << " --------------------------------------------- " << endl;
	LocalCheckingStats tmp_local_stats(_max_num_prime_factors);
	_local_stats_vec.resize(_local_cover.multiplicity(), tmp_local_stats);
	cout << " Finished initializing the local cover statistics data for " << _max_num_prime_factors << " local covers " << endl << endl << endl;

	  
	// Make the quaternary theta series 
	cout << " Computing the quaternary theta series: " << endl;
	cout << " ---------------------------------------------------- " << endl;
	_GetQuaternaryTheta(qf_data);
	cout << " Finished computing the quaternary theta series " << endl << endl << endl;
	

	
	// Find the squarefree exceptions
	cout << " Checking the eligible squarefree numbers: " << endl;
	cout << " ----------------------------------------- " << endl;
	_CheckSquarefree_New(qf_data);
	cout << " Finished checking the eligible squarefree numbers " << endl << endl << endl;
	
	
	// Find all exceptions
	cout << " Checking the eligible square factors: " << endl;
	cout << " ------------------------------------- " << endl;
	_CheckSquareFactors();
	cout << " Finished checking the eligible square factors " << endl << endl << endl;
      }
      
    }
    else
      // If the cusp constant is zero, then we still should write that there are no (square-free) exceptions!
      _WriteSquarefree(qf_data);
    
  }  
  
} 




////////////////////////
// Secret Constructor
////////////////////////
void Representability::_MakeWithUniformCover(const QF_Datafiles & qf_data, const Matrix_mpz & QQ, const double & cusp_const) {

  // TO DO:
  //   Check that Q is 4x4, symmetric, and positive definite

  // Set the Overflow Bound
  const double OVERFLOW_BOUND = 2000;
  overflow = false;    // Default value needed for when the cusp_const is zero.

  // Store the (pre-computed) cuspidal upper bound
  _passed_cusp_const = cusp_const;
  cout << " Passed in the cusp constant " << cusp_const << endl << endl;

  
  // Compute the level and the character of QQ
  _BigForm = QQ;
  _level = QQ.QFLevel().get_si();
  _character = SquarefreePart(QQ.Determinant()).get_si();
  
  // TO LOOK AT:
  //   unsigned long Precision = (unsigned long) ceil(Theta_Precisions[k]);


  // Try to read the (square-free) exceptions from a file first  <--- **** NEED TO CHANGE THIS WEHN WE CHECK SQUARE FACTORS TOO ***
  if (_ReadSquarefree(qf_data) == true)
    _used_squarefree_exception_file == true;
  else {
    _used_squarefree_exception_file == true;
    
    
    // Skip all of the work if there are no cusp forms ! =)
    if (_passed_cusp_const != 0) {
      
      // Find the F4 Bound for QQ
      cout << " Finding the F4 Bound: " << endl;
      cout << " --------------------- " << endl;
      _GetF4Bound(qf_data);
      cout << " Finished finding the F4 Bound " << endl << endl << endl;
      
      // Check for an overflow... (when F4_bound > OVERFLOW_BOUND)
      if (_F4_bound > OVERFLOW_BOUND) {
	overflow = true;
	cout << " Representability Overflow: The bound " << _F4_bound << " is > " << OVERFLOW_BOUND << endl;
      }
      else {
	overflow = false;
	
	
	// Find the eligible primes 
	cout << " Finding the eligible primes: " << endl;
	cout << " ---------------------------- " << endl;
	_GetEligiblePrimes(qf_data);
	cout << " Finished finding the eligible primes " << endl << endl << endl;
	
	// Check that there are actually numbers to check!
	if (_max_num_prime_factors != 0) {
	  cout << " _max_num_prime_factors = " << _max_num_prime_factors << endl << endl;
	  
	  // Construct a uniform local cover for QQ
	  cout << " Making the uniform local cover: " << endl;
	  cout << " ------------------------------- " << endl;
	  _local_cover.Make("uniform", QQ, qf_data.Local_Cover_Dir);
	  cout << " Finished making the uniform local cover " << endl << endl << endl;
	  
	  
	  // Determine the exception set for our cover
	  cout << " Making the set of exceptions for the uniform local cover: " << endl;
	  cout << " --------------------------------------------------------- " << endl;
	  set<mpz_class> Exceptions_from_Uniform_Cover;
	  /*  FIX THIS LATER...
	  Exceptions_from_Uniform_Cover = _Make_Uniform_Cover_Exceptions(_local_cover, 0, qf_data);
	  */
	  mpz_class largest_exception = 0;
	  if (Exceptions_from_Uniform_Cover.empty() == false)
	    largest_exception = *(Exceptions_from_Uniform_Cover.rbegin());
	  cout << " Finished making the set of exceptions for the uniform local cover " << endl;
	  cout << " The largest exception is " << largest_exception << endl;
	  cout << " There are " << Exceptions_from_Uniform_Cover.size() << " exceptions." << endl << endl << endl;
	  
	  
	  // Check all exceptions to see if they're represented by Q               // <--- *** SHOULD BE REPLACED BY AN EXCEPTIONS ROUTINE! *** 
	  cout << " Finding the squarefree exceptions of Q: " << endl;
	  cout << " --------------------------------------- " << endl;
	  PowerSeries<mpz_class> big_theta = _BigForm.ComputeTheta(largest_exception.get_ui() + 2);
	  set<mpz_class> Exceptions;
	  for (set<mpz_class>::iterator i = Exceptions_from_Uniform_Cover.begin(); i != Exceptions_from_Uniform_Cover.end(); i++)
	    if (big_theta[(*i).get_ui()] == 0) 
	      Exceptions.insert(*i);
	  cout << " Finished finding the squarefree exceptions of Q. " << endl;	
	  cout << " They are: " << Exceptions << endl << endl << endl;
	  
	  // Copy the set to the squarefree exceptions vector
	  for (set<mpz_class>::iterator i = Exceptions.begin(); i != Exceptions.end(); i++)
	    _squarefree_exceptions.push_back(*i);
	  
	  
	  
	  // Deal with square factors (i.e. the numbers m*p^2) for p|m and p<=11
	  
	  
	  
	  // Deal with anisptropic primes
	  
	  
	  
	  
	  // Copy the set to the all_exceptions vector
	  for (set<mpz_class>::iterator i = Exceptions.begin(); i != Exceptions.end(); i++)
	    _all_exceptions.push_back(*i);	
	  
	  
	}	
	else {
	  // Here we still need to check the number 1, and some squares...         ***** TO DO *****

	  // If there are numbers to check, then we still should write that there are no (square-free) exceptions!
	  _WriteSquarefree(qf_data);
	}

	
      }
	
    }
    else
      // If the cusp constant is zero, then we still should write that there are no (square-free) exceptions!
      _WriteSquarefree(qf_data);
    
  }  
  
} 



///////////////////////////////////////////////////////////////////////////////////////////////
// Find the time (in some units) needed to compute the set of exceptions for the minimal cover 
///////////////////////////////////////////////////////////////////////////////////////////////
double Representability::_thetatime_for_minimal_local_cover(const LocalSplitCoverings & cover) const {

  double temp;  // Temporary variable to accumulate the time estimates

  for(unsigned long i=0; i < cover.multiplicity(); i++)
    for(unsigned long j=0; j < cover.size(i); j++)
      temp += _ThetaPrecision(cover.Get_Value(i,j)).get_d();    // Note: _ThetaPrecision() assumes 5 attempts!
  
  return temp;

}


///////////////////////////////////////////////////////////////////////////////////////////////
// Find the time (in some units) needed to compute the set of exceptions for the uniform cover, assuming the largest exception is E.
///////////////////////////////////////////////////////////////////////////////////////////////
double Representability::_thetatime_for_uniform_local_cover(const LocalSplitCoverings & cover, const mpz_class & E) const {

  // Sanity check that the cover is actually uniform
  for(unsigned long i=0; i < cover.multiplicity(); i++)
    for(unsigned long j=0; j < cover.size(i); j++)
      assert(cover.Get_Value(i,j) == cover.Get_Value(0,0));

  // Declare some convenient variables
  double temp;                                               // Temporary variable to accumulate the time estimates
  const double M = _Upper_bound_for_eligible_numbers();
  const double d = cover.Get_Value(0,0);

  // Compute the precision needed for the theta computations    (See 1/15/05 notes for details.) 
  const double B = (pow(d * E.get_d(), 0.6666) * pow(M, 0.3333) + 0.25 * pow(d * E.get_d(), 1.3333) * pow(M, -0.3333)); 

  // Compute the number of forms in the local cover
  long number_of_forms = 0;
  for(unsigned long i=0; i < cover.multiplicity(); i++)
    number_of_forms += cover.size(i);
  
  // Return the time estimate    (See 1/15/05 notes for details.) 
  return number_of_forms * B;

}



////////////////////////////////////////////////////////////////////////////////////////////////
// Find the maximum size for the maximal exception E needed for the uniform theta computation to be faster than the minimal one.
// (Based on the assumption that M >> E.  See 1/15/05 notes for details.) 
////////////////////////////////////////////////////////////////////////////////////////////////
mpz_class Representability::_uniform_maximum_exception_target(const LocalSplitCoverings & cover) const {

  // Sanity check that the cover is actually uniform
  for(unsigned long i=0; i < cover.multiplicity(); i++)
    for(unsigned long j=0; j < cover.size(i); j++)
      assert(cover.Get_Value(i,j) == cover.Get_Value(0,0));

  // Declare some useful constants
  const double T_old = _thetatime_for_minimal_local_cover(cover);
  const double M = _Upper_bound_for_eligible_numbers();
  const double d = cover.Get_Value(0,0);

  // Compute the number of forms in the local cover
  long number_of_forms = 0;
  for(unsigned long i=0; i < cover.multiplicity(); i++)
    number_of_forms += cover.size(i);
  

  // DIAGNOSTIC 
  cout << "  --> Computed the exception target to be: " << pow(T_old / number_of_forms, 1.5) / (d * sqrt(M)) << endl;
  cout << "  -->            As an mpz_class, this is: " << mpz_class(pow(T_old / number_of_forms, 1.5) / (d * sqrt(M))) << endl;



  return mpz_class(pow(T_old / number_of_forms, 1.5) / (d * sqrt(M)));

}



///////////////////////////////////////////////////////////////////////
// Chooses whether to use the minimal cover or uniform cover method...
//////////////////////////////////////////////////////////////////////
string Representability::_choose_checking_method(const QF_Datafiles & qf_data, const double & F4) const {

  // Check if the bound is large enough to consider the uniform method
  if (F4 <= 500) 
    return "minimal";

  // Make a minimal split local cover, and get a theta time estimate.
  LocalSplitCoverings minimal_cover("minimal", _BigForm);
  const double Minimal_time = _thetatime_for_minimal_local_cover(minimal_cover);
  
  // Make a uniform split local cover, and get an estimate for the largest allowed exception.
  LocalSplitCoverings uniform_cover("uniform", _BigForm);
  const mpz_class Exception_Max = _uniform_maximum_exception_target(uniform_cover);


  // See if this holds for our uniform cover.
  set<mpz_class> Exceptions_from_Uniform_Cover;
  /*  FIX THIS LATER...
  Exceptions_from_Uniform_Cover = _Make_Uniform_Cover_Exceptions(uniform_cover, Exception_Max, qf_data);
  */



  //                                                   ********* FINISH THIS ROUTINE! **********

}


/*


///////////////////////////////////////////////////////////////////////////////////////
// Find the set of exceptions so long as the largest is less than the desired bound.
///////////////////////////////////////////////////////////////////////////////////////

set<mpz_class> Representability::_Make_Uniform_Cover_Exceptions(const LocalSplitCoverings & cover, 
								const mpz_class & desired_exception_bound,
								const QF_Datafiles & qf_data) const {
								

  // Sanity check that the cover is actually uniform
  for(unsigned long i=0; i < cover.multiplicity(); i++)
    for(unsigned long j=0; j < cover.size(i); j++)
      assert(cover.Get_Value(i,j) == cover.Get_Value(0,0));

  // Declare some useful constants
  const double M = _Upper_bound_for_eligible_numbers();                   // Perhaps we should pass this???
  const double d = cover.Get_Value(0,0);


  // Make a vector to store the cumulative list of exceptions...
  vector<mpz_class> cumulative_bound_list;                // Stores the bound used for each iteration.
  vector< set<mpz_class> > cumulative_exceptions_list;    // Stores the cleaned exceptions list for each iteration of the bound.
  set<mpz_class> temporary_exceptions_list;               // Stores the exceptions-in-progress for each iteration.

  // Run through each form of the cover, until we hit our running_bound 
  // or the largest exception hits the desired_exception_bound.
  mpz_class starting_bound = 0;
  mpz_class running_bound = mpz_class(2 * sqrt(M*d) + d);
  mpz_class list_size = 0;                       // Says the size of the cumulative list
  mpz_class biggest_exception = 0;
  mpz_class A;
  bool done_flag = false;                        // Says when we have finished checkign the exceptions! =)


  // Compute the number of forms in the local cover
  unsigned long number_of_forms = 0;
  for(unsigned long i=0; i < cover.multiplicity(); i++)
    number_of_forms += cover.size(i);
  


  while (((desired_exception_bound == 0) || (biggest_exception <= desired_exception_bound)) && (done_flag == false)) {
    
    // Find the (incremental) exceptions of each ternary up to running_bound
    vector< set<mpz_class> > ternary_excep_vector(number_of_forms);
    for(unsigned long i=0; i < cover.multiplicity(); i++)
      for(unsigned long j=0; j < cover.size(i); j++) {
	ternary_excep_vector[i] = TernaryExceptions(cover.Get_Form(i,j), running_bound, qf_data.Ternary_Exceptions_Dir).GetExceptions(); 
	ternary_excep_vector[i].erase(ternary_excep_vector[i].begin(), ternary_excep_vector[i].upper_bound(starting_bound));     // Takes the subset >= starting_bound...
      }


    
    
    
    // Make a combined list of possible exceptions for the cover
    for(unsigned long i=0; i < number_of_forms; i++)
      temporary_exceptions_list.insert(ternary_excep_vector[i].begin(), ternary_excep_vector[i].end());

    
    // Check that each possible exception is either an exception for, or not locally represented by, every ternary.
    set<mpz_class> cover_exceptions;
    for(set<mpz_class>::iterator i = temporary_exceptions_list.begin(); i != temporary_exceptions_list.end(); i++) {
      bool exception_flag = true;
      
      // Check if the number is represented by some (other) ternary
      for(unsigned long j=0; j < number_of_forms; j++) 
	if((cover.Get_Local_Conditions(j).IsLocallyRepresented(*i) == true) && (ternary_excep_vector[j].find(*i) == ternary_excep_vector[j].end()))
	  exception_flag = false;
      
      // Keep the actual exceptions for the cover
      if (exception_flag == true)
	cover_exceptions.insert(*i);
           
    }


    // Add the exceptions/bound to the cumulative list of cover exceptions/bounds
    cumulative_bound_list.push_back(running_bound);
    cumulative_exceptions_list.push_back(cover_exceptions);
    list_size++;


    // DIAGNOSTIC
    cout << endl;
    cout << " ternary_excep_vector.size() = " << ternary_excep_vector.size() << endl;
    cout << " ternary_excep_vector = " << ternary_excep_vector << endl;
    cout << " starting_bound = " << starting_bound << endl;
    cout << " running_bound = " << running_bound << endl;
    cout << " temporary_exceptions_list = " << temporary_exceptions_list << endl;
    cout << " cover_exceptions = " << cover_exceptions << endl;
    cout << endl;



    // Update the largest exception
    if (cover_exceptions.empty() == false)
      biggest_exception = *(cover_exceptions.rbegin());


    // Find the smallest A which proves the theorem, assuming that we have the overall largest exception.
    // (Based on our writeup...)                                                                        ***** CHECK THIS!!! *****
    A = 1;
    double AA = A.get_d();
    while ( (8*AA*(AA+1)*(AA+1)) * sqrt(M/d) < (biggest_exception + 1) * (biggest_exception + 1) ) {      // Is this ok??

      // DIAGNOSTIC
      cout << " Looking for the minimal A value: " << endl;
      cout << " E = " << biggest_exception << endl;
      cout << " (E+1)^2 = " << ((biggest_exception + 1) * (biggest_exception + 1)) << endl;
      cout << " 8*A*(A+1)^2 * sqrt(M/d) = " << ((8*AA*(AA+1)*(AA+1)) * sqrt(M/d)) << endl;
      cout << " (8*AA*(AA+1)*(AA+1)) * sqrt(M/d) < (biggest_exception + 1) * (biggest_exception + 1) = " 
	   << ((8*AA*(AA+1)*(AA+1)) * sqrt(M/d) < (biggest_exception + 1) * (biggest_exception + 1)) << endl;
      cout << endl;

      AA++;

    }

    // Update the starting and running_bound for the next iteration, based on the new A
    mpz_class new_running_bound = mpz_class(2 * (A+1) * sqrt(M*d) + d * (A*A - 1));
    if (new_running_bound <= running_bound)
      done_flag = true;
    else {
      starting_bound = running_bound;
      running_bound = new_running_bound;
    }
    
  }
  
  
  // Return {0} if we fail... =(
  if ((desired_exception_bound != 0) && (biggest_exception > desired_exception_bound)) {
    cout << endl;
    cout << " Failed to find exceptions small enough... =( " << endl;
    cout << "   Desired Exception Bound: " << desired_exception_bound << endl;
    cout << "   Biggest Exception found: " << biggest_exception << endl;
    cout << endl;

    set<mpz_class> fail_set;
    fail_set.insert(0);
    return fail_set;
  }


  // Otherwise, return the (conglomerated) exception set
  set<mpz_class> good_set;
  for(unsigned long i=0; i < cumulative_bound_list.size(); i++) 
    good_set.insert(cumulative_exceptions_list[i].begin(), cumulative_exceptions_list[i].end());

  cout << endl;
  cout << " Success! =) " << endl;
  cout << "   Using A = " << A << " and bound = " << running_bound << "." << endl;
  cout << "   The Largest Exception is: " << biggest_exception << endl;
  cout << "   The cover has " << good_set.size() << " exceptions." << endl;
  cout << "   The exception set is " << good_set << endl;
  cout << endl;

  return good_set;
}

*/






////////////////
// Destructor
///////////////
Representability::~Representability() {

  // Do Nothing... =)

}



///////////////////////
// Get the Exceptions 
///////////////////////
vector<mpz_class> Representability::GetExceptions() const {

  return _all_exceptions;

}


//////////////////////////////////////////////////
// (Somehow) Get the F4_bound for the form...
///////////////////////////////////////////////////
void Representability::_GetF4Bound(const QF_Datafiles & qf_data) {

  // Get the Eisenstein coefficient lower bound
  mpq_class Eis_lower_bound;
  cout << " Computing the theoretical Eisenstein coefficient bound... " << endl;
  Eis_lower_bound = _BigForm.GetEisensteinLowerBound(qf_data.Eis_Lower_Bound_Dir); 
  cout << " Finished computing the theoretical Eisenstein coefficient bound... " << endl << endl;




  // Sanity Check (part 1): Compute the against the first few Eisenstein Coefficients
  mpq_class Eis_numerical_bound;  
  cout << " Computing the (approximate) numerical Eisenstein coefficient bound... " << endl;
  Eis_numerical_bound = _BigForm.NumericalEisensteinLowerBound(qf_data.Eis_Dir);
  cout << " Finished computing the (approximate) numerical Eisenstein coefficient bound " << endl << endl;



  //  /*
  // Summary
  cout << "   The theoretical lower bound is: " << Eis_lower_bound << " which is about " << (Eis_lower_bound.get_d() * 1.0) << endl;
  cout << "   The numerical lower bound is: " << Eis_numerical_bound << " which is about " << (Eis_numerical_bound.get_d() * 1.0) << endl;
  cout << endl;
  //  */




  // Sanity Check (part 2): Compare the bound againt the first few Eisenstein coefficients
  if (Eis_lower_bound <= Eis_numerical_bound)
    cout << " The theoretical lower bound is less than the numerical lower bound from the coefficients! =)" << endl;
  else {
    cout << " ERROR: There's a problem since the theoretical lower bound is larger than what we observe in the coefficients! =( " << endl;
    assert(0 == 1);
  }
  cout << endl;


  // Get the cuspidal upper bound
  double Cusp_upper_bound;
  if (_passed_cusp_const < 0) {
    Cusp_upper_bound = 100;              // ********** THIS NEEDS TO BE CHANGED TO COMPUTE THE CORRECT CONSTANT!!! *************
  } else
    Cusp_upper_bound = _passed_cusp_const;


  // Make the F4 upper bound
  _F4_bound = Cusp_upper_bound / Eis_lower_bound.get_d();
  cout << " The F4 upper bound is: " << _F4_bound << endl;


  // Read this from a data file, or compute it by hand...
 
} 


//////////////////////////////////////////////////////////
// Get the list of possible primes with F4(p) <= Bound
// (Based on PrimeListFromBoundOptimized().)
/////////////////////////////////////////////////////////
void Representability::_GetQuaternaryTheta(const QF_Datafiles & qf_data) {


  // Initialize the quaternary theta function to the correct precision
  PowerSeries<mpz_class> tmp_series(_BIGFORM_THETA_PRECISION);
  _big_theta = tmp_series;



  // Get the Quaternary Theta series:
  // --------------------------------

  // If we use persistent data then read or write the primes and bounds, otherwise just compute them.
  if (true) {                                                //   <----- TO DO: Can replace this later by: (qf_data.Keep_Quaternary_Theta == true) 
      
    // Create the output filename: "quaternary_theta__FORM__PRECISION.txt"
    string quaternary_theta_filename;
    quaternary_theta_filename = qf_data.Theta_Dir + "quaternary_theta__" + _BigForm.QF_String() + "__" + MakeString(_BIGFORM_THETA_PRECISION) +".txt";  
        
    /*
    // DIAGNOSTIC
    cout << "\n Filename: " << quaternary_theta_filename.str() << "X" << endl << endl;  
    */

    // Read the theta directory if it exists
    bool read_flag = _ReadQuaternaryTheta(quaternary_theta_filename);
    
    // Otherwise, make the lists and write them to files
    if (read_flag == false) {
      /*
      _big_theta = _BigForm.GetMagmaThetaSeries("", _BIGFORM_THETA_PRECISION);   // Use MAGMA Routine
      */
      _big_theta = Theta_PARI_1_new(_BigForm, _BIGFORM_THETA_PRECISION);         // Use my (PARI-based) routine!
      _big_theta.WriteSeries(quaternary_theta_filename.c_str()); 
    }

  }
  else 
      _big_theta = _BigForm.GetMagmaThetaSeries("", _BIGFORM_THETA_PRECISION);

}




////////////////////////////////////////////////////////////////////////////
// Check if the prime and F4(p) files exist, and if so then read them in.
////////////////////////////////////////////////////////////////////////////
bool Representability::_ReadQuaternaryTheta(const string & quaternary_theta_filename) {

  // Check if the file already exists...
  bool File_exists = FileExists(quaternary_theta_filename);

  // If it exists, then read it in and exit. =)
  if (File_exists == true) {
    _big_theta.ReadSeries(quaternary_theta_filename.c_str());

    cout << " Read the quaternary theta function from the datafile." << endl;

    // If they have the same number of elements then exit, otherwise proceed below...
    if (_big_theta.Precision() == _BIGFORM_THETA_PRECISION)       // WARNING: These tests may be moot since the precision doesn't come from the file...
      return true;
    else {
      cout << " _ReadQuaternaryTheta() Error:  The precision of the theta function is not 10,000... =( ... so we'll have to regenerate it..." << endl;
      return false;
    }
  }
  
  // otherwise, continue and regenerate the files...
  else {
    cout << " _ReadQuaternaryTheta() Error:  The precision of the theta function is not 10,000... =( ... so we'll have to regenerate it..." << endl;
    return false;
  }


  // Sanity Check
  cout << " _ReadQuaternaryTheta() Error: We shouldn't be here!  Abort!" << endl;
  abort();

}







//////////////////////////////////////////////////////////
// Get the list of possible primes with F4(p) <= Bound
// (Based on PrimeListFromBoundOptimized().)
/////////////////////////////////////////////////////////
void Representability::_GetEligiblePrimes(const QF_Datafiles & qf_data) {


  /*
  // DIAGNOSTIC
  cout << " Entering PrimeListFromBound with: " << endl;
  cout << "       bound = " << _F4_bound << endl;
  cout << "       level = " << _level << endl;
  cout << "     chi_top = " << _character << endl;
  */


  // Get the Eligible Primes:
  // ------------------------

  // If we use persistent data then read or write the primes and bounds, otherwise just compute them.
  if (qf_data.Keep_Eligible_Primes == true) {
      
    // Create the output filenames: 
    //    primes:  "eligible_primes__bound_N_chitop.data2"
    //     F4(p):  "eligible_F4(p)__bound_N_chitop.data2"
    string primefilename, F4filename;
    primefilename = qf_data.Eligible_Prime_Dir + "eligible_primes__b=" + MakeString(_F4_bound) + 
      "_N=" + MakeString(_level) + "_chi=" + MakeString(_character) + "__.data2"; 
    F4filename = qf_data.Eligible_Prime_Dir + "eligible_F4(p)__b=" + MakeString(_F4_bound) + 
      "_N=" + MakeString(_level) + "_chi=" + MakeString(_character) + "__.data2"; 
        
    /*
    // DIAGNOSTIC
    cout << "\n Filename: " << primefilename.str() << "X" << endl << endl;  
    cout << "\n Filename: " << F4filename.str() << "X" << endl << endl;  
    */

    // Read eligible primes and their F4(p) constants if they exist
    bool read_flag = _ReadEligiblePrimes(primefilename, F4filename);
    
    // Otherwise, make the lists and write them to files
    if (read_flag == false) {
      _ComputeEligiblePrimes();   
      _WriteEligiblePrimes(primefilename, F4filename);
    }

  }
  else 
    _ComputeEligiblePrimes();   



  // Find the maximum number of prime factors allowed in an eligible square-free number
  _FindMaxNumPrimeFactors();


  // Print the results:
  cout << " There are " << _plist.size() << " eligible prime numbers." << endl;
  long VVnum = 5;
  if (_plist.size() < 5) 
    VVnum = _plist.size();
  if (_plist.size() > 0) {
    cout << " The last " << VVnum << " eligible primes and their F4(p) values are: " << endl;
    PrintTailV(_plist, VVnum);
    PrintTailV(_f4list, VVnum);    
  }

}



////////////////////////////////////////////////////////////////////////////
// Check if the prime and F4(p) files exist, and if so then read them in.
////////////////////////////////////////////////////////////////////////////

bool Representability::_ReadEligiblePrimes(const string & primefilename, const string & F4filename) {

  // Check if these files already exist...
  bool Files_exist = (FileExists(primefilename) && FileExists(F4filename));

  // If they exist, then read them in and exit. =)
  if (Files_exist == true) {
    _plist = ReadVector_long(primefilename.c_str());
    _f4list = ReadVector_double(F4filename.c_str());

    cout << " Read the eligible primes and the F4 list from the datafiles." << endl;

    // If they have the same number of elements then exit, otherwise proceed below...
    if (_plist.size() == _f4list.size())
      return true;
    else {
      cout << " _ReadEligiblePrimes() Error:  The # of primes and the # of F4(p) don't match, so we'll have to regenerate them..." << endl;
      return false;
    }
  }
  
  // otherwise, continue and regenerate the files...
  else {
    cout << " _ReadEligiblePrimes() Error:  Can't open eligible prime and F4(p) files, so we'll have to regenerate them... " << endl;
    return false;
  }


  // Sanity Check
  cout << " _ReadEligiblePrimes() Error: We shouldn't be here!  Abort!" << endl;
  abort();

}



////////////////////////////////////////////////////////////////////////////////////
// Write the results of the eligible primes and their F4 Constants F4(p) to files
////////////////////////////////////////////////////////////////////////////////////

bool Representability::_WriteEligiblePrimes(const string & primefilename, const string & F4filename) const {


  // Open the output files and check they're opened correctly
  ofstream primefile, F4file; 
  primefile.open(primefilename.c_str());
  F4file.open(F4filename.c_str());

  if (! primefile.is_open())
    { cout << "_WriteEligiblePrimes(): Error opening output file " << primefilename; exit (1); }
  if (! F4file.is_open())
    { cout << "_WriteEligiblePrimes(): Error opening output file " << F4filename; exit (1); }
  

  // Write the primes and their F4(p) values
  for(long i=0; i < _plist.size(); i++) {
    primefile << _plist[i] << ", ";
    F4file << _f4list[i] << ", ";
  }


  // Remove the trailing ", ", insert an EOF, and close the files
  long file_ptr;
  file_ptr = primefile.tellp();
  primefile.seekp(file_ptr - 2);
  primefile.put(' ');
  primefile.put(EOF);
  primefile.close();

  file_ptr = F4file.tellp();
  F4file.seekp(file_ptr - 2);
  F4file.put(' ');
  F4file.put(EOF);
  F4file.close();

}






///////////////////////////////////////////////////////////////
// Define the function F4primeiso which computes F4(p)
//   (assuming p is prime and isotropic)
///////////////////////////////////////////////////////////////
double Representability::_F4PrimeIsoExact(const long p) const {

  double f4;
  double pp = double(p);  // Warning: This conversion is necessary to prevent a subtle error...(see below)
  f4 = sqrt(pp) * 0.5;
  if ((KroneckerSymbol(_character, p) == -1) && ((_level % p) != 0))
    f4 = f4 * (pp - 1) / (pp + 1);
  //    f4 = f4 * (pp*pp - pp + 1) / (pp*pp + 1);   // Why was this here?? =| (J.H. 7/12/05)


  return f4;
}




//////////////////////////////////////////////////////////
// Compute the list of possible primes with F4(p) <= Bound,
// and sort them by the size of their F4(p).
// (Based on PrimeListFromBoundOptimized().)
/////////////////////////////////////////////////////////

void Representability::_ComputeEligiblePrimes() {

  // Sanity check that there are at least 100 primes is our list
  assert(100 < Big_Prime_List.size() - 1);


  // Increase the F4_bound to account for all F4(p) < 1 in our search for eligible primes! =)
  double Modified_F4_bound_for_primes = _F4_bound;  
  vector<unsigned long> small_primes;
  small_primes.push_back(2); small_primes.push_back(3);
  small_primes.push_back(5); small_primes.push_back(7);

  for (int i=0; i < small_primes.size(); i++) {

    unsigned long p = small_primes[i];
    double F4_val = _F4PrimeIsoExact(p);
    if (F4_val < 1) 
      Modified_F4_bound_for_primes = Modified_F4_bound_for_primes / F4_val;

  }


  // /*
  // DIAGNOSTIC:
  cout << " Computing eligible primes: " << endl;
  cout << "    Original F4 bound for eligible numbers= " << _F4_bound << endl;  
  cout << "    Modified F4 bound for eligible primes= " << Modified_F4_bound_for_primes << endl;  
  PrintTime();
  // */



  // Set some local variables to look for eligible primes
  unsigned long Fast_F4_transition = 10000;    // Compute the F_4(p) bound exactly for all primes p < 10,000
  unsigned long prime_ptr = 0;
  unsigned long p = Big_Prime_List[prime_ptr];      // This sets p=2
  double pp, F4_val, F4_val_min=0;


  // Find all eligible primes
  while ( (F4_val_min <= Modified_F4_bound_for_primes) && (prime_ptr < Big_Prime_List.size() - 1) ) {
        
    // Compute the (approximate) lower bound for the next F4(p)
    pp = double(p);                        // Warning: This conversion is necessary to prevent a subtle overflow error...
    F4_val_min = sqrt(pp) * 0.5;
    F4_val_min = F4_val_min - F4_val_min/pp;  // This is the approximate F4 value which is below F4(p). =)
    
    // Compute the next F4(p) value (depending on the size of p)
    if (p <= Fast_F4_transition) 
      F4_val = _F4PrimeIsoExact(p);
    else
      F4_val = F4_val_min;

    // If p is eligible, add p and F4(p) to the list
    if (F4_val <= Modified_F4_bound_for_primes) {
      _plist.push_back(p);
      _f4list.push_back(F4_val);
    }

    // Increment the prime
    prime_ptr = prime_ptr + 1;
    p = Big_Prime_List[prime_ptr];
    
  }


  // ERROR TRAP:
  if (prime_ptr >= Big_Prime_List.size() - 1) {
    cout << "ERROR in _ComputeEligiblePrimes(): We ran out of primes... =(" << endl;
    cout << "    Trying to use a larger list of primes." << endl;

    // Try to use the next prime file and repeat the computation
    Use_Next_Primefile();
    _ComputeEligiblePrimes();
  }                                                                                                        


  // DIAGNOSTIC
  /*
  cout << " The _plist is : " << _plist << endl;
  cout << " The _f4list is : " << _f4list << endl;
  */
  cout << endl;
  cout << " Finished making the eligible primes, now we're sorting them! " << endl;
  PrintTime();

  
  // Check if there's more than one eligible prime to sort
  if  (_plist.size() > 1) {

    // Sort the eligible primes in order of increasing F4(p)
    bool sorted_flag = false;
    while (sorted_flag == false) {
      sorted_flag = true;
      unsigned long tmp_p; double tmp_f4;
      cout << " Making a (bubble sort) pass at ordering the primes by F4(p)..."; 
      PrintTime();
      for(long i=0; i < min(10000, _plist.size() - 1); i++)
	if (_f4list[i+1] < _f4list[i]) {

	  // OUTPUT
	  cout << " Swapping p=" << _plist[i] << " with q=" << _plist[i+1]
	       << ":    F4(" << _plist[i] << ")=" << _f4list[i] << " and F4(" << _plist[i+1] << ")=" << _f4list[i+1] << endl;

	  /*
	  // DIAGNOSTIC
	  cout << " previous sorted_flag = " << sorted_flag << endl;
	  cout << " i = " << i << endl;
	  cout << " p[i] = " << _plist[i] << "   F4[i] = " << _f4list[i] << endl;
	  cout << " p[i+1] = " << _plist[i+1] << "   F4[i+1] = " << _f4list[i+1] << endl;
	  cout << " (_f4list[i+i] < _f4list[i]) = " << (_f4list[i+i] < _f4list[i]) << endl;
	  cout << endl;
	  */
	  
	  // Swap primes
	  tmp_p = _plist[i];
	  _plist[i] = _plist[i+1];
	  _plist[i+1] = tmp_p;
	  
	  // Swap values
	  tmp_f4 = _f4list[i];		 
	  _f4list[i] = _f4list[i+1];		 		 
	  _f4list[i+1] = tmp_f4;		 		 

	  // Reset the flag
	  sorted_flag = false;
	}
    }

  }

  // DIAGNOSTIC
  cout << " Finished sorting the eligible primes! " << endl;
  PrintTime();



  // OUTPUT
  cout << " Exiting _ComputeEligiblePrimes() with: " << endl;
  cout << "      Number of Eligible primes: " << _plist.size() << endl;
  long output_tail_size = 5;
  if (_plist.size() >= output_tail_size) {
    cout << "      The first "<< output_tail_size << " pairs (p, F4(p)) are: " << endl << "        ";
    for (long ii=0; ii < output_tail_size; ii++)
      cout << " (" << _plist[ii] << ", " << _f4list[ii] << ") ";
    cout << endl;
    
    cout << "      The last "<< output_tail_size << " pairs (p, F4(p)) (in reverse order) are: " << endl << "        ";
    for (long ii=_plist.size() - 1; ii > (_plist.size() - output_tail_size - 1); --ii)
      cout << " (" << _plist[ii] << ", " << _f4list[ii] << ") ";
    cout << endl;
  }
  PrintTime();


  /*
  // DIAGNOSTIC:
  cout << "           p = " << p << endl;
  cout << "       F4(p) = " << F4_val << endl;
  cout << "   prime_ptr = " << prime_ptr << endl;
  */  

}




// Find out the maximal number of prime factors allowed in an eligible squarefree number 
// (Note: This is complicated since we need to keep track of 
// the order of the primes, but we'd really like to keep
// track of the order of the F4 bounds...  
void Representability::_FindMaxNumPrimeFactors() {

  long num = 0;
  double temp_F = 1;

  //  TO DO: THIS CAN BE CLEANED UP A LITTLE SINCE THE PRIMES ARE IN INCREASING F4 ORDER!


  // First take the product of all the F4(p) < 1 to find the smallest it can be
  // (This can only happen if p = 2, 3, 5, or 7!)
  for(unsigned long i=0; i<=3; i++)
    if ((i < _f4list.size()) && (_f4list[i] <= 1)) {
      temp_F = temp_F * _f4list[i];  
      num++;
      cout << " The product after " << num << " factors is: " << temp_F << endl;
    }

  // If this is still too big, then there are no eligible numbers! =)
  if (temp_F > _F4_bound)
    _max_num_prime_factors = 0;    

  else {  
    
    // Keep adding more primes until we're too big...
    unsigned long i=0;
    while ((temp_F <= _F4_bound) && (i < _f4list.size())) {  
      if (_f4list[i] > 1) {
	temp_F = temp_F * _f4list[i];
	num++;
	cout << " The product after " << num << " factors is: " << temp_F << endl;
      }

      i++;
    }  

    num--;    // Decrement if we overflow!
    
    cout <<  "We can have square-free numbers with at most " << num << " prime factors." << endl;
    
    // Set the class variable
    _max_num_prime_factors = (long) num;
    
  }
}






/*!  \brief Computes the precision for the ternary theta series 
 * needs to check our quaternary form.
 *
 * This computes an upper bound for the precision needed for the 
 * ternary form to check that our given 4 variable form is represented.  
 * This bound only depends on the splitting coefficient \f$d\f$ in the 
 * decomposition \f$Q = T \oplus d w^2\f$.
 *
 * The proof of our upper bound for \f$T\f$ is
 * \f[
 * B > F_4(T) = \prod_p F_4(p) 
 *            = \prod_p \left( \frac{\sqrt{p}}2 \cdot adj(p) \right)
 * \hfill \Rightarrow \hfill
 *  T = \prod_p p \leq \left( \frac{2^{\# p|T} B} {adj(p)} \right)^2. 
 * \f] 
 *
 * To bound the maximal difference, we assume that \f$ T = dx^2 \f$
 * (which is the worst case) and that we will attempt to find the
 * minimal difference \f$t\f$ times.  Then 
 * \f[
 * T - d(x-t)^2 = dx^2 - d(x-t)^2 = d(2tx - t^2) < 2tdx = 2t\sqrt{Td}, 
 * \f]
 * since \f$ x = \sqrt{\frac Td} \f$.
 *
 *
 *
 *
 * \todo This should take an input number \f$t\f$ for the number of tries.
 *
 */


//////////////////////////////////////////////////////////////////////////////////
// Computes an upper bound for the size of the largest eligible squarefree number
//    Note: This also sets the # of prime factors...
/////////////////////////////////////////////////////////////////////////////////

double Representability::_Upper_bound_for_eligible_numbers() const {

  // Sanity Check that _max_num_prime_factors has been set
  // assert( _max_num_prime_factors > 0);                            // NOTE: This may be zero for some forms... so we removed it...


  // TO DO:
  //   Get rid of the prime list here...


  //  /*
  // DIAGNOSTIC
  cout << " Entering _Upper_bound_for_eligible_numbers: " << endl;
  cout << " Using B = " << _F4_bound 
       << "  N = " << _level 
       << "  chi_top = " << _character
    //       << "  diag_coeff = " << diag_coeff   
       << "  max_num_prime_factors = " << _max_num_prime_factors << endl;  
  //  */



  // Deal with the bounds of zero gracefully
  if (_F4_bound == 0)
    return 0;


  // Make a small list of primes
  const long prime_list[] = 
    {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 
     101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 
     197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 
     311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 
     431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 
     557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 
     661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 
     809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 
     937, 941, 947, 953, 967, 971, 977, 983, 991, 997}; 
  const long Prime_List_Length = 168;

  

  // The global list of primes we can use. =)
  // (Declared in main.cc)
  //extern const vector<long> Big_Prime_List;  




  // Find the product of the first N adjustment factors (N = Max allowed number of primes) 
  long ct = 0;
  double adjustment = 1;
  long i = 0; long p = Big_Prime_List[i];    
  if (_character != 1) {        // Only check for an improvement if the character is non-trivial! =)  This is ok because _character is made square-free! 
    while ((ct < _max_num_prime_factors) && (i < Big_Prime_List.size())) {        // Use < or <= ???
      if (((_level % p) != 0) && (KroneckerSymbol(_character, p) == -1)) {
	/*
	// Describes the primes where we make a slight improvement...
	cout << "Making an adjustment for p = " << p << endl;
	*/
	adjustment = adjustment * (double(p-1) / double(p+1));
	ct++;
      }
      
      i++; p = Big_Prime_List[i];      
    }
  }




  // *** QUESTION: What's the deal with the character being 1 here???

  // Error Checking:  Check we haven't used all of the primes or chi_top is 1. =)
  if ((i == Big_Prime_List.size()) && (_character != 1)){
    cout << "ERROR2 in ThetaPrecision:  We've used up all of the primes < 1000. " << endl;
    cout << "  Your bound must have been ridiculously huge, or the primes are biased! " << endl;
    exit(1);
  }


  // Compute the needed precision 
  double T_bound;
  T_bound = pow(2.0, 1.0 * _max_num_prime_factors) * _F4_bound / adjustment;
  T_bound = T_bound * T_bound;

  /*
    cout << " 1/adjustment is " << (1/adjustment) << endl;
    cout << " T_bound is " << T_bound << endl;
  */

  cout << " The largest eligible number must be less than " << T_bound << endl <<endl;

  return T_bound;

}



////////////////////////////////////////////////////////////////////
// Compute an upper bound for the precision needed in the ternary //
// theta series to check the eligible square-free numbers.        //
////////////////////////////////////////////////////////////////////

mpz_class Representability::_ThetaPrecision(const long diag_coeff, const long iteration_depth) const {

  cout << endl << endl << " Computing the Ternary Precision " << endl;

  // Sanity Checks
  assert(iteration_depth >= 1);

  // Compute the precision for depth 1
  double Diff_bound;
  Diff_bound = 2 * floor(sqrt(_Upper_bound_for_eligible_numbers())) * diag_coeff;

  /*
    cout << " Diff_bound is " << Diff_bound << endl;
  */

  // INSERT A FACTOR OF 5 TO GET THE CORRECT PRECISION FOR A DEPTH 5 HUNT! -- See 8/4/04 notes.
  Diff_bound = DEFAULT_CHECKING_DEPTH * Diff_bound; 
  cout << " Starting with a default number-checking depth of " << DEFAULT_CHECKING_DEPTH 
       << " which gives the ternary precision " << Diff_bound << endl;


  // We will start with depth 5, and then double the precision for each iteration
  for(long i=1; i<=iteration_depth-1; i++) {
    cout << " Doubling the ternary precision from " << Diff_bound << " to " << (2 * Diff_bound) << endl; 
    Diff_bound *= 2;
  }

  cout << endl << " Finished computing the Ternary Precision " << endl << endl;

  return mpz_class( ceil(Diff_bound) );

}




/////////////////////////////////////////////////////////////////////////////////////////////
// Make the vector of boolean ternary theta functions for the ternaries in our local cover   // <----  THIS ROUTINE SHOULD BE DELETED!!!!
/////////////////////////////////////////////////////////////////////////////////////////////

bool Representability::_MakeTernaryThetas(const QF_Datafiles & qf_data, const long iteration_depth) {
  
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


  // Clear any old boolean ternary theta functions
  _ternary_theta.clear();


  // Loop through each local cover
  for(long i=0; i < _local_cover.multiplicity(); i++) {

    vector<boolean_theta> tmp_vec;   // Temporarily store the thetas for each local cover


    // Loop through each form in this local cover
    for(long j=0; j < _local_cover.size(i); j++) {
      
      // Compute the precision needed for each form, and get its boolean theta function
      mpz_class Precision = _ThetaPrecision(_local_cover.Get_Value(i,j), iteration_depth);  
      
      // Check that the precision is withing acceptable bounds  (Precision <= 1 million)
      if (Precision > 1000000)
	return false;   
      
      
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
      
      
      // Read the theta function, and if that fails then compute and write it!
      if ( temp_theta.read(theta_filename.c_str()) == false ) {      
	cout << " Can't read the boolean theta function, so we're recomputing it!" << endl;
	
	
	// Use our approximate (PARI-based) routine to make the boolean theta function directly
	temp_theta = Theta_PARI_2_new_Approximate(_local_cover.Get_Form(i,j), Precision);


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
      tmp_vec.push_back(temp_theta);
      
    }
      
    // Add the thetas for this cover to the big list of boolean ternary theta functions
    _ternary_theta.push_back(tmp_vec);

  }
  
  
  // (Frivolous) Sanity check to ensure we have as many boolean theta functions as forms.
  assert (_ternary_theta.size() == _local_cover.multiplicity());
  for(unsigned long i=0; i < _local_cover.multiplicity(); i++)
    assert (_ternary_theta[i].size() == _local_cover.size(i));
  
  
  // This means the ternary precision is acceptable! =)
  return true;
  
}





/////////////////////////////////////////////////////////////////////////////////////////////
// Make the vector of boolean ternary theta functions for the ternaries in our local cover
/////////////////////////////////////////////////////////////////////////////////////////////

void Representability::_UpdateTernaryThetas(const QF_Datafiles & qf_data) {
  
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


  // Some sanity checks
  assert(_local_cover.multiplicity() >= _ternary_theta.size());    // Check there are more local covers than computed theta functions.
  for(long i=0; i < _ternary_theta.size(); i++) {             // Check that the number of computed theta functions agrees 
    cout << " i = " << i << endl;
    cout << " _local_cover.size(i) = " << _local_cover.size(i) << endl;
    cout << " _ternary_theta[i].size() = " << _ternary_theta[i].size() << endl;
    assert(_local_cover.size(i) == _ternary_theta[i].size());      //   with the number of forms in the local covers.
  }

  //  /*
  // DIAGNOSTIC:
  cout << endl << " Entering _UpdateTernaryThetas: " << endl;
  cout << "   _ternary_theta.size() = " << _ternary_theta.size() << endl;
  cout << "   _local_cover.multiplicity() = " << _local_cover.multiplicity() << endl;
  cout << endl;
  //  */


  // Loop through each new local cover
  for(long i=_ternary_theta.size(); i < _local_cover.multiplicity(); i++) {

    vector<boolean_theta> tmp_vec;   // Temporarily store the thetas for each local cover


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


	// Use our approximate (PARI-based) routine to make the boolean theta function directly
	temp_theta = Theta_PARI_2_new_Approximate(_local_cover.Get_Form(i,j), Precision);

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
      tmp_vec.push_back(temp_theta);
      
    }
      
    // Add the thetas for this cover to the big list of boolean ternary theta functions
    _ternary_theta.push_back(tmp_vec);

  }
  
  
  // (Frivolous) Sanity check to ensure we have as many boolean theta functions as forms.
  assert (_ternary_theta.size() == _local_cover.multiplicity());
  for(unsigned long i=0; i < _local_cover.multiplicity(); i++)
    assert (_ternary_theta[i].size() == _local_cover.size(i));
  
  
}





/////////////////////////////////////////////////////////////
// Write the (square-free) exceptions if we're supposed to.
/////////////////////////////////////////////////////////////

void Representability::_WriteSquarefree(const QF_Datafiles & qf_data) const {

  if (qf_data.Keep_Exceptions == true) {

    // Write the squarefree exceptions:
    // --------------------------------
    
    // Make the filename
    string excep_filename;
    excep_filename = qf_data.Exceptions_Dir + "squarefree__" + _BigForm.QF_String() + ".txt";

    // Write the (square-free) exceptions
    ofstream excep_stream;
    excep_stream.open(excep_filename.c_str());
    FileOutput(_squarefree_exceptions, excep_stream);
    excep_stream.close();


    // Write the squarefree overflows:
    // -------------------------------
        
    // Make the filename
    string overflow_filename;
    overflow_filename = qf_data.Overflows_Dir + "squarefree_overflow__" + _BigForm.QF_String() + ".txt";

    // Write the (square-free) exceptions
    ofstream overflow_stream;
    overflow_stream.open(overflow_filename.c_str());
    FileOutput(_squarefree_overflows, overflow_stream);
    overflow_stream.close();

  }

}



/////////////////////////////////////////////////////////////
// Read the (square-free) exceptions if we can.
/////////////////////////////////////////////////////////////

bool Representability::_ReadSquarefree(const QF_Datafiles & qf_data) {
  
  // Read the squarefree exceptions:
  // -------------------------------
    
  // Make the filename
  string excep_filename;
  excep_filename = qf_data.Exceptions_Dir + "squarefree__" + _BigForm.QF_String() + ".txt";
  
  // Check if we can read the exceptions
  ifstream excep_stream;
  excep_stream.open(excep_filename.c_str());
  if (excep_stream.fail())
    return false;  // Read failed... =(

  /*
  // DIAGNSOTIC
  cout << " First the square_free exceptions are " << _squarefree_exceptions << endl;
  */  


  // Read the (square-free) exceptions, if we can.
  FileInput(_squarefree_exceptions, excep_stream);
  excep_stream.close();


  /*
  // DIAGNSOTIC
  cout << " Now the square_free exceptions are " << _squarefree_exceptions << endl;
  */



  // Read the squarefree overflows:
  // ------------------------------
  
  // Make the filename
  string overflow_filename;
  overflow_filename = qf_data.Overflows_Dir + "squarefree_overflow__" + _BigForm.QF_String() + ".txt";
  
  // Write the (square-free) exceptions
  ifstream overflow_stream;
  overflow_stream.open(overflow_filename.c_str());
  if (overflow_stream.fail())
    return false;  // Read failed... =(

  // Read the (square-free) overflows, if we can.
  FileInput(_squarefree_overflows, overflow_stream);
  overflow_stream.close();


  return true;     // Read succeeded! =)
  
}




////////////////////////////////////////////////////////////////////////////////////////////
// Incrementally generates eligible squarefree numbers, and checks their representability.
/////////////////////////////////////////////////////////////////////////////////////////////

void Representability::_CheckSquarefree(const QF_Datafiles & qf_data) {


  // Make a counter to count squarefree numbers
  // and lists of exceptions and overflows...
  mpz_class Counter;
  Counter = 0;
  vector<mpz_class> exception_list, overflow_list, overflow_differences, overflow_depths;

  vector<mpz_class> T_list;

  int Max_num = _max_num_prime_factors;


  long CHECKING_DEPTH = 100;    // This is the deepest we will check for an exception!
                                // Is we check this far, and 


  // Check if the big form is locally universal (once and for all)
  bool Bigform_universal_flag = _local_cover.Big_Local().IsUniversal();



  // Make the ternary theta function for the maximum allowed depth
  cout << "   There are " << _local_cover.multiplicity()  << " split local coverings." << endl;
  cout << "   Computing the relevant boolean ternary theta series for pre-existing local covers: " << endl;
  _UpdateTernaryThetas(qf_data);
  cout << "   Finished computing the relevant boolean ternary theta series for pre-existing local covers! " << endl << endl << endl;
  




  

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

	// Check these squarefree numbers for execptions
	_NewestCheckNumbers_by_number(qf_data, T_list, CHECKING_DEPTH);	


	T_list.resize(0);

	cout << "    Cumulative results: " << endl;
	cout << "      Checked " << Counter << " numbers " << endl;
	cout << "      Number of exceptions: " << _squarefree_exceptions.size() << endl;
	cout << "      Number of overflows: " << _squarefree_overflows.size() << endl;
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
  sort(_squarefree_exceptions.begin(), _squarefree_exceptions.end());  


  // Write the (square-free) exceptions if we're supposed to.
  _WriteSquarefree(qf_data);


  // ------ Find the square-free exceptions of the original form (not just the sublattice) --------

  // Compute the theta function of Q up to the largest exception
  //   This is the last exception:  exception_list[exception_list.size()-1]


  
  // Write the exceptions (for both) to files
  

}






// Check the representability of a vector of eligible numbers.
void Representability::_NewestCheckNumbers_Subroutine(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const int Max_Depth, 
						      vector<mpz_class> & possible_exception_list) {

  // Clear the vector to return
  possible_exception_list.clear();


  // Check each of the numbers in num
  for (unsigned long i=0; i< num.size(); i++) {

    bool done_flag = false;                      // This says when to finish using ternaries to check representability of a number
    bool repn_flag = false;                      // This says when a number is represented by some split local covering form
    unsigned long cover_number = 0;
    unsigned long ternary_form_number = 0;
    const long MAXIMUM_NUMBER_OF_LOCAL_COVERS = 10;          // This is the maximum of local covers we'll use before giving up! =) 
    const long BIGGEST_EXCEPTION = 50000;                    // This is the bound for how large an exception must be to look for another local cover!
    
    // If we run out of local split covers, then generate more up to some preset maximum, when we give up!
    while ((done_flag == false) && (repn_flag == false) && (_ternary_theta.size() <= MAXIMUM_NUMBER_OF_LOCAL_COVERS)) {

      
      // Until the number is represented or we're out of forms, run through all forms
      while ((repn_flag == false) && (cover_number < _ternary_theta.size())) {
	
	bool overflow_flag = false;  // This detects an overflow of the difference m - dx^2 for the ternary form

	while ((repn_flag == false) && (overflow_flag == false) && (ternary_form_number < _ternary_theta[cover_number].size())) {


	  long depth = 0;
	  
	  // Then run through all allowed depths, checking the differences for representability
	  while ((repn_flag == false) && (depth < Max_Depth)) {
	    
	    // Find the w to use when subtracting d*(w)^2 
	    mpz_class w;
	    long d = _local_cover.Get_Value(cover_number, ternary_form_number);
	    w = mpz_class(floor(sqrt(num[i].get_d() / d))) - depth;    // Check this is correct...
	    if (w < 0)  
	      w = 0;	
	    
	    // Find the difference
	    mpz_class diff;
	    diff = num[i] - mpz_class(d) * w * w;
	    
	    /*
	    // DIAGNOSTIC
	    if (num[i] == 290) {
	    cout << " ternary form number = " << ternary_form_number << endl;
	    cout << " ternary_theta.size() = " << _ternary_theta.size() << endl;
	    cout << " depth = " << depth << endl;
	    cout << "   d = " << d << endl;
	    cout << "   i = " << i 
	    << "   num[i] = " << num[i] << endl;
	    cout << "   w = " << w 
	    << "   diff =  " << diff << endl;
	    cout << endl;
	    }
	    */
	    
	    
	    // Sanity Checks
	    if (diff < 0) {
	      cout << "Error in computing the difference for num[" << i <<"] = "<< num[i] << endl;
	      exit(1);
	    }
	    if (diff.fits_uint_p() == false) {
	      cout << "The difference " << diff << " for num = " << num[i] << " is too large to handle!" << endl;
	      exit(1);
	    }
	    
	    
	    
	    // Check if the difference is represented by the current ternary:
	    // --------------------------------------------------------------
	    
	    // Check that each difference is in range
	    if (diff <= _ternary_theta[cover_number][ternary_form_number].precision()) {
	      
	      // If so, check if the ternary theta series represents the difference
	      if (_ternary_theta[cover_number][ternary_form_number].get_value(diff.get_ui()) == true)    
		repn_flag = true;


	      
	      // cout << " ------ Looking at the number " << num[i] << endl;
	      
	      /*
	      // DIAGNOSTIC
	      if (num[i] == 290) {
	      cout << " --> repn_flag = " << repn_flag << endl;
	      
	      cout << endl;
	      for (unsigned long mm=0; mm<=30; mm++)
	      cout << mm << " -> " << _ternary_theta[ternary_form_number].get_value(mm) << endl;
	      cout  << endl;
	      }
	      */
	      
	    }
	    else 
	      overflow_flag = true;   // Say we overflow the precision of this form...
	    

	    // Finished checking with this depth, so increment
	    depth++;
	  }
	  
	  // Finished checking with this ternary form, so increment
	  ternary_form_number++;
	}
	
	// Finished checking with this local cover, so increment
	cover_number++;
	ternary_form_number = 0;
      }

      
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
      
 
    } 
    

    
    /*
    // DIAGNOSTIC
    if (num[i] == 290) {
      cout << " repn_flag = " << repn_flag << endl;
      cout << " Is it locally represented? " << _local_cover.Big_Local().IsLocallyRepresented(num[i]) << endl;
      cout << endl;
    }
    */
    
    
    
    // If the number is not represented using any ternary (and not an overflow), 
    // but it is locally represented, then add it to the exception list. =)
    if ((repn_flag == false) && (_local_cover.Big_Local().IsLocallyRepresented(num[i]) == true))
      possible_exception_list.push_back(num[i]);	      
    
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
void Representability::_NewestCheckNumbers_by_number(const QF_Datafiles & qf_data, const vector<mpz_class> & num, const int Max_Depth) {


  // Here is where we deal with overflows and decide which possible 
  // square-free exceptions and overflows are genuine exceptions!
  // --------------------------------------------------------------


  // Declare some constants and local variables
  vector<mpz_class> possible_exception_list, overflow_list; 
  const long BIGGEST_EXCEPTION = 10000;
  const long EXCEPTION_LIMIT = 1000000;              // Don't continue if there are exceptions this big at the end.
  const long BIGGEST_NUMBER_OF_EXCEPTIONS = 100;
  vector<mpz_class> new_eligible_list;
  mpz_class max_eligible;
  long iter_depth = 1;
  bool ternary_precision_ok = true;  // Says when the ternary precision has gotten too large!

  // Check the vector of numbers
  _NewestCheckNumbers_Subroutine(qf_data, num, Max_Depth, possible_exception_list);




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











void Representability::_CheckSquareFactors() {

  // DIAGNSOTIC
  cout << endl;
  cout << " ---> The square-free exceptions are: " << _squarefree_exceptions << endl;
  cout << endl;


  // TO CHANGE: Temporarily ignore the square factors...
  _all_exceptions = _squarefree_exceptions;

}










