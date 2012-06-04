

// Empty Constructor 
LocalSplitCoverings::LocalSplitCoverings() {

  // Do Nothing... This is needed for some future class declarations...

}


// Constructor which finds the local split covering for QQ 
LocalSplitCoverings::LocalSplitCoverings(const string & cover_type, const Matrix_mpz & QQ) {

  Make(cover_type, QQ, "");

}

// Constructor which finds the local split covering for QQ 
LocalSplitCoverings::LocalSplitCoverings(const string & cover_type, const Matrix_mpz & QQ, const string & LocalCoverDir) {

  Make(cover_type, QQ, LocalCoverDir);

}



// Finds the local split covering for QQ 
void LocalSplitCoverings::Make(const string & cover_type, const Matrix_mpz & QQ, const string & LocalCoverDir) {
  
  // Sanity Checks:
  assert((cover_type == "minimal") || (cover_type == "uniform"));    // Check that the cover types are valid.
    // Q is non-zero, symmetric, and non-singular



  // Set the big form, find its local conditions, and set the cover_type
  _big_form = QQ;
  _big_local.Make(QQ);
  _cover_type = cover_type;


  // Make the filename if we pass a non-empty directory
  if (LocalCoverDir != "") 
    _filename_string = LocalCoverDir + "Local_Cover__" + QQ.QF_String() + "__" + _cover_type + ".txt";
  else 
    _filename_string.erase();  


  // DIAGNOSTIC:
  cout << " LocalCoverDir = " << LocalCoverDir << endl;
  cout << " _filename_string = " << _filename_string << endl;

  
  // Get the Local Cover:       // Note: Here the empty directory "" is interpereted as a directive not to save the data!
  // --------------------
  __Largest_size_used = 0;     // Set that we haven't computed any theta_vectors yet! =)
  if ((_filename_string.empty() == false) && (FileExists(_filename_string) == true))    // Check if we've already computed a local cover
    _ReadLocalSplitCoverings(); 
  else {

    // Compute the split local cover of the appropriate type (and write it if an output_string is set) 
    if (_cover_type == "minimal") {
      _IncrementMinimalLocalSplitCoverings();          // Compute a minimal local split covering, then write it to a file (if we set a filename)
      
    }
    else if (_cover_type == "uniform") {
      cout << "LocalSplitCoverings::Make() ERROR: There is currently no support for uniform coverings... =(" << endl;
      exit(1);
      //_FindUniformLocalSplitCoverings();          // Compute a uniform local split covering 
    }
    else {
      cout << " Error in LocalSplitCoverings::Make:  The cover type " << cover_type << " is not allowed.";
      assert(0==1);
    }
    
  } 
  

  // DIAGNOSTIC
  cout << " This is the local cover: " << endl;
  cout << " ------------------------ " << endl;
  cout << " _big_form: " << endl << _big_form << endl;
  cout << " _big_local: " << endl << _big_local << endl;
  cout << " _form_list: " << endl << _form_list << endl;
  cout << " _d_list: " << endl << _d_list << endl;
  cout << " _local_list: " << endl << _local_list << endl;
  cout << " _cover_type: " << endl << _cover_type << endl;
  cout << endl << endl << endl;  
  



  // Sanity Check:
  // Check that the lengths of all the vectors are the same...
  assert(_d_list.size() == _form_list.size());
  assert(_local_list.size() == _form_list.size());

}


// Add an additional local cover to a pre-defined local covering
void LocalSplitCoverings::IncrementLocalSplitCoverings() {

  // Some basic sanity checks to ensure the cover has already been defined
  assert((_cover_type == "minimal") || (_cover_type == "uniform"));
  assert(_big_form.NumRows() == 4);
  assert(_big_form.IsSymmetric() == true);

  // Call the correct increment routine
  if (_cover_type == "minimal") {
    _IncrementMinimalLocalSplitCoverings();          // Compute a minimal local split covering       
  }
  else if (_cover_type == "uniform") {
    cout << "LocalSplitCoverings::IncrementLocalSplitCoverings() ERROR: There is currently no support for uniform coverings... =(" << endl;
    exit(1);
    //_FindUniformLocalSplitCoverings();          // Compute a uniform local split covering 
  }

}




// Empty Destructor
LocalSplitCoverings::~LocalSplitCoverings(){
  // Do nothing...
}



/////////////////////////////////////////
// Gives the number of local coverings //
/////////////////////////////////////////
unsigned long LocalSplitCoverings::multiplicity() const{
  return _form_list.size();  
};


////////////////////////////////////////////////////////////////////////
// Gives the number of forms in the i-th local covering (where i>=0). //
////////////////////////////////////////////////////////////////////////
unsigned long LocalSplitCoverings::size(const long i) const{
  return _form_list[i].size();  
};




////////////////////////////////////////////////////////////////////////////////////////
// Gives the "big" form, whose local conditions are covered by codimension 1 subforms
////////////////////////////////////////////////////////////////////////////////////////
Matrix_mpz LocalSplitCoverings::Big_Form() const{
  return _big_form;  
};



////////////////////////////////////////////////////////////////////////////////////////
// Gives the local conditions of "big" form
////////////////////////////////////////////////////////////////////////////////////////
LocalConditions LocalSplitCoverings::Big_Local() const{
  return _big_local;  
};



////////////////////////////////////////////////////////
// Gives the i-th subform in the split local covering
////////////////////////////////////////////////////////
Matrix_mpz LocalSplitCoverings::Get_Form(const long cover_number, const long form_number) const{

  if ((cover_number >= 0) && (cover_number < _form_list.size()) 
      && (form_number >= 0) && (form_number < _form_list[cover_number].size()))
    return _form_list[cover_number][form_number];
  else {
    cout << "Error in LocalSplitCoverings::Get_Form(): Trying to access an out of range element..." << endl;
    exit(1);
  }
  
};


/////////////////////////////////////////////////////////////////////////////////////////////////////
// Gives the i-th value (of the vector perpendicular to the i-th form) in the split local covering
/////////////////////////////////////////////////////////////////////////////////////////////////////
long LocalSplitCoverings::Get_Value(const long cover_number, const long form_number) const{

  if ((cover_number >= 0) && (cover_number < _d_list.size()) 
      && (form_number >= 0) && (form_number < _d_list[cover_number].size()))
    return _d_list[cover_number][form_number];
  else {
    cout << "Error in LocalSplitCoverings::Get_Value(): Trying to access an out of range element..." << endl;
    assert(0==1);
  }
  
};


/////////////////////////////////////////////////////////////////////////////////////////////////////
// Gives the i-th value (of the vector perpendicular to the i-th form) in the split local covering
/////////////////////////////////////////////////////////////////////////////////////////////////////
LocalConditions LocalSplitCoverings::Get_Local_Conditions(const long cover_number, const long form_number) const{

  if ((cover_number >= 0) && (cover_number < _local_list.size()) 
      && (form_number >= 0) && (form_number < _local_list[cover_number].size()))
    return _local_list[cover_number][form_number];
  else {
    cout << "Error in LocalSplitCoverings::Get_Local_Conditions(): Trying to access an out of range element..." << endl;
    assert(0==1);
  }
  
};






///////////////////////////////////////////////
// Writes a LocalSplitCoverings to the ostream
///////////////////////////////////////////////
void LocalSplitCoverings::_FileOutput(ostream & out) const {

  // Write the opening braces " <<< "
  out << " <<< " << endl;

  // Write the cover_type
  out << _cover_type;
  out << endl << endl;

  // Write the big form
  FileOutput(_big_form, out); 
  out << endl << endl;

  // Write the local conditions for the big form
  FileOutput(_big_local, out);   
  out << endl << endl;

  // Write the (n-1)-dim'l forms in the local covering 
  FileOutput(_form_list, out);   
  out << endl << endl;

  // Write the local conditions for the (n-1)-dim'l forms in the local covering 
  FileOutput(_local_list, out);   
  out << endl << endl;

  // Write the 1-dim'l forms in the local covering 
  FileOutput(_d_list, out);   
  out << endl;

  // Write the list of used indices for these local covers
  FileOutput(__used_indices, out);   
  out << endl;

  // Write the closing braces " >>> "
  out << " >>> " << endl;

}



///////////////////////////////////////////////
// Reads a LocalSplitCoverings from the istream
///////////////////////////////////////////////
void LocalSplitCoverings::_FileInput(istream & in) {

  // Declare some temporary variables
  Matrix_mpz tmp_big_form;
  LocalConditions tmp_big_local;
  vector< vector<Matrix_mpz> > tmp_form_list;
  vector< vector<long> > tmp_d_list;
  vector< vector<LocalConditions> > tmp_local_list;
  set< vector<long> > tmp__used_indices;
  char ch;

  // Read the opening braces " <<< "
  for (long i=1; i<=3; i++) {
    in >> ch;
    assert(ch = '<');
  }

  // Read the cover_type
  in >> _cover_type;
  assert((_cover_type == "minimal") || (_cover_type == "uniform"));    // Check that the cover types are valid.

  // Read the big form
  FileInput(tmp_big_form, in); 

  // Read the local conditions for the big form
  FileInput(tmp_big_local, in);   

  // Write the (n-1)-dim'l forms in the local covering 
  FileInput(tmp_form_list, in);   

  // Write the local conditions for the (n-1)-dim'l forms in the local covering 
  FileInput(tmp_local_list, in);   

  // Write the 1-dim'l forms in the local covering 
  FileInput(tmp_d_list, in);   

  // Write the list of used indices for these local covers
  FileInput(tmp__used_indices, in);   


  // Read the closing braces " >>> "
  for (long i=1; i<=3; i++) {
    in >> ch;
    assert(ch = '>');
  }


  // Some Sanity checks
  assert( tmp_d_list.size() == tmp_form_list.size() );
  assert( tmp_d_list.size() == tmp_local_list.size() );


  // Copy the temporary variables
  _big_form = tmp_big_form;
  _big_local = tmp_big_local;
  _form_list = tmp_form_list;
  _local_list = tmp_local_list;
  _d_list = tmp_d_list;
  __used_indices = tmp__used_indices;
}













/*
// TODO: Make this active...

////////////////////////////////////////////////////
// This is the overloaded (global) version for <<
////////////////////////////////////////////////////
ostream & operator<<(ostream & out, const LocalSplitCoverings & local) {    
  return local.Print(out);
}
*/









/////////////////////////////////////////////////////////////////////////
// Finds all small vectors for QQ, ordered by the numbers they represent.
/////////////////////////////////////////////////////////////////////////
vector< vector< vector<long> > > LocalSplitCoverings::_ThetaVectors(const Matrix_mpz & QQ, const unsigned long bound) const {

  unsigned long Theta_Precision = bound;

  
  // 1.Compute the theta function to precision 70, and then keep track of
  // each of the vectors in a vector of vectors.

  /*
  // Print the ternary form and its level
  cout << endl;
  cout << " The form is: " << endl; 
  cout << QQ << endl;
  cout << " The level of the form is " << QQ.QFLevel() << endl;
  */


  const long n = QQ.NumRows();
  
  double Cholesky_temp[n * (n+1)];    // This should get declared in CholeskyDecomposition(...)!
  CholeskyDecomposition(QQ, Cholesky_temp);

  /*  
  cout << " Computing the theta function " << endl;
  PrintTime();

  cout << endl << "Entering FastBinaryTheta" << endl;
  cout << " Using precision = " << precision() << endl;
  */

  // ERROR CHECKING: Check that C+1 fits in an unsigned long.


  // Make the vector of vectors which have a given value
  // (So theta_vec[i] will have all vectors v with Q(v) = i.)
  vector< vector<long> > empty_vec_list;
  vector< vector< vector<long> > > theta_vec;
  theta_vec.resize(Theta_Precision + 1, empty_vec_list);




  
  // Initialize Q with zeros
  double Q[n+1][n+1];
  for (long i=1; i<=n; i++)
    for (long j=1; j<=n; j++) 
      Q[i][j] = 0;


  // Copy the Cholesky array into Q
  long cholesky_counter = 0;
  for (long i=1; i<=n; i++)
    for (long j=i; j<=n; j++) {
      Q[i][j] = Cholesky_temp[cholesky_counter];
      cholesky_counter++;
    }
 


  /*
  // Print the floating-point matrix Q
  cout << "Using Q = " << endl;
  for (long i=1; i<=n; i++) {
    for (long j=1; j<=n; j++)
      cout << Q[i][j] << " ";
    cout << endl;
  }
  */


  // 1. Initialize
  long i = n;
  vector<double> T(n+1, 0);  // Note: We use n+1 so we can index the entries as 1 --> n
  vector<double> U(n+1, 0);
  T[i] = (double) Theta_Precision;
  U[i] = 0;
  double Z;
  vector<long> L(n+1, 0);
  vector<long> x(n+1, 0);


  // 2. Compute bounds
  Z = sqrt(T[i] / Q[i][i]);
  L[i] = long(floor(Z - U[i]));  // Check this is ok...
  /*
  cout << " L[i] float gives :     " << floor(Z - U[i].get_d()) << endl;
  cout << " L[i] mpz_class gives : " << L[i] << endl;
  */

  x[i] = long(ceil(-Z - U[i]) - 1);  // Check this is ok...
  /*
  cout << " x[i] float gives :     " << ceil(-Z - U[i].get_d() -1) << endl;
  cout << " x[i] mpz_class gives : " << x[i] << endl;
  cout << endl;
  */
  
  
  bool done_flag = false;
  double Q_val_double;
  unsigned long Q_val;                 // WARNING: Still need a good way of checking overflow for this value...

  
  // Big loop which runs through all vectors
  while (done_flag == false) {

    // Loop through until we get to i=1 (so we defined a vector x)
    do {
      
      // 3a. Main loop
      x[i] = x[i] + 1;
      while (x[i] > L[i]) {
	i = i + 1;
	x[i] = x[i] + 1;
      }
      
      // 3b. Main loop
      if (i>1) {
	/*
	cout << " i = " << i << endl;
	cout << " T[i] = " << T[i] << endl;
	cout << " Q[i][i] = " << Q[i][i] << endl;
	cout << " x[i] = " << x[i] << endl;
	cout << " U[i] = " << U[i] << endl;
	cout << " x[i] + U[i] = " << (x[i] + U[i]) << endl;	
	cout << " T[i-1] = " << T[i-1] << endl;		
	*/
	T[i-1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i]);
	/*
	cout << " T[i-1] = " << T[i-1] << endl;		
	cout << endl;
	*/
	i = i - 1;
	U[i] = 0;
	for(long j=i+1; j<=n; j++)
	  U[i] = U[i] + Q[i][j] * x[j];
	
	// Now go back and compute the bounds...
	// 2. Compute bounds
	Z = sqrt(T[i] / Q[i][i]);
	L[i] = long(floor(Z - U[i]));  
	x[i] = long(ceil(-Z - U[i]) - 1);  
      }
      
    } while (i > 1);
    
    
    // 4. Solution found (This happens when i=1)
    /*
    cout << " x = " << x << endl;
    cout << " Q_val = Q(x) = " << Q_val << endl;
    */
    Q_val_double = Theta_Precision - T[1] + Q[1][1] * (x[1] + U[1]) * (x[1] + U[1]);
    Q_val = (unsigned long) round(Q_val_double);

    //    cout << " Float = " << Q_val_double << "   Long = " << Q_val << "  XX " << endl;
    /*
    cout << " The float value is " << Q_val_double << endl;
    cout << " The associated long value is " << Q_val << endl;
    */

    if (Q_val <= Theta_Precision) {
      vector<long> y;
      //      cout << " Have vector " << x << " with value " << Q_val;
      y = x;
      y.erase(y.begin());               // Trim the irrelevant 0th entry, 
      theta_vec[Q_val].push_back(y);    // then add the vector.
      //      cout << " which we replace by " << y << endl;
      //      cout << endl;
    }


 
    // 5. Check if x = 0, for exit condition. =)
    long j=1;
    done_flag = true;
    while (j<=n) {
      if (x[j] != 0)
	done_flag = false;
      j++;
    }    
  }


  
  cout << " Leaving ThetaVectors()" << endl;
  

  return theta_vec;

}











//////////////////////////////////////////////////////////////////////////////
// Finds the (ternary) form orthogonal to the vector vec inside of BigForm.
//////////////////////////////////////////////////////////////////////////////
Matrix_mpz LocalSplitCoverings::_FindComplementarySublattice(const Matrix_mpz & QQ, const vector<long> & vec) const {

  assert(QQ.NumRows() == QQ.NumCols());

  /*
  // DIAGNOSTIC
  cout << " Entering FindComplementarySublattice() with " << endl
       << "   vec = " << vec << endl
       << "   and Bigform = " << endl << BigForm;

  cout << " The value of BigForm[v] is " << BigForm.EvaluateQuadratic(vec) << endl << endl;
  */

  // Make a new Matrix_mpz Q for QQ
  Matrix_mpz Q;
  Q = QQ;
  long n = Q.NumRows();

  // Find the first non-zero component of vec, and call it nz  (Note: 0 <= nz < n)
  long nz = 0;
  while ((nz < n) && (vec[nz] == 0))
    nz++;

  // Abort if vec is the zero vector
  if (nz == n) {
    cout << "Error in FindComplementarySublattice(): vec is the zero vector!" << endl;
    exit(1);
  }

  // Make the change of basis matrix
  Matrix_mpz new_basis;  // This is the change of basis matrix
  new_basis = IdentityMatrix(n);
  for(long i=1; i<=n; i++)   // Put vec in the nz-th column
    new_basis(i, nz+1) = vec[i-1];

  /*
  cout << "The new basis matrix is: " << endl
       << new_basis << endl;
  new_basis.Transpose();
  cout << "The new basis matrix transposed is: " << endl 
       << new_basis << endl;
  cout << "The new basis matrix transposed (again) is: " << endl 
       << new_basis.GetTranspose() << endl;
  */

  // Change Q (to Q1) to have vec as its nz-th basis vector 
  Matrix_mpz Q1(n,n);
  /*
  cout << " Q = " << endl << Q; 
  cout << " Q1 = " << endl << Q1; 
  cout << " new_basis = " << endl << new_basis; 
  cout << " new_basis.GetTranspose() = " << endl << new_basis.GetTranspose(); 
  */
  Q1 = new_basis.GetTranspose() * Q * new_basis;
  /*
  cout << " Q = " << endl << Q; 
  cout << " Q1 = " << endl << Q1; 
  cout << " new_basis = " << endl << new_basis; 
  cout << " new_basis.GetTranspose() = " << endl << new_basis.GetTranspose(); 
  cout << "\n\n";
  */

  // Pick out the value Q(vec) of the vector
  mpz_class d = Q1(nz+1, nz+1);


  // For each row/column, perform elementary operations to cancel them out.
  for(long i=1; i <= n; i++){

    // Don't change the pivot row/column!
    if (i != nz+1) {

      // Check if the (i,nz+1)-entry is divisible by d,
      // and stretch its row/column if not.
      if (Q1(i,nz+1) % d != 0) 
	Q1 = MultiplySymmetric(Q1, d / GCD(d, Q1(i,nz+1)), i);
      
      // Now perform the (symmetric) elementary operations to cancel out the (i,1) entries/
      Q1 = AddSymmetric(Q1, -Q1(i,nz+1) / GCD(d, Q1(i,nz+1)), i, nz+1);

    }
  }

  // Check that we're done!
  bool done_flag = true;
  for(long i=1; i <= n; i++) {
    if((i != nz+1) && (Q1(nz+1,i) != 0))
      done_flag = false;
  }
  if (done_flag == false) {
    cout << "\n Error, Dude! =0 " << endl;

    // DIAGNOSTIC
    cout << " Entering FindComplementarySublattice() with " << endl
	 << "   vec = " << vec << endl
	 << "   and QQ = " << endl << QQ;
    
    cout << " The value of BigForm[v] is " << QQ.EvaluateQuadratic(vec) 
	 << " while the value of d is " << d << endl;
    
    cout << " The change of basis matrix is: " << endl << new_basis;
    cout << " And Q (should be the same as BigForm) is: " << endl << new_basis;

    
    cout << "ERROR IN FindDecomposableTernary.....of1 -- The form is not orthogonal... =(" << endl;
    cout << " Q1 is " << endl << Q1;
    cout << " and nz = " << nz << endl; 
  }

  
  // Return the complementary matrix
  Matrix_mpz T;
  T.SetDims(n-1,n-1);
  for(long i=1; i <= n; i++)
    for(long j=1; j <= n; j++)
      if ((i != nz+1) && (j != nz+1))
      
	if ((i < nz+1) && (j < nz+1))
	  T(i,j) = Q1(i,j);
	else if ((i < nz+1) && (j > nz+1))
	  T(i,j-1) = Q1(i,j);
	else if ((i > nz+1) && (j < nz+1))
	  T(i-1,j) = Q1(i,j);
	else if ((i > nz+1) && (j > nz+1))
	  T(i-1,j-1) = Q1(i,j);
  

  return T;

}






// Want to keep track of the local conditions for the splitting
// sublattices, and order them from best (close to universal) to
// worst. Then choose a minimal set of lattices which gives the local
// conditions of the original form.

////////////////////////////////////////////////////////////////////////////////////
// Finds a local split covering of _big_form, perpendicular to vectors of minimal length.
////////////////////////////////////////////////////////////////////////////////////
void LocalSplitCoverings::_IncrementMinimalLocalSplitCoverings() {

  /*
  // Sanity check that the big form is locally universal
  assert(_big_local.IsUniversal() == true);
  */

  // Set some preferences
  const long BATCH_INCREMENT = 5;    // Find all vectors of length 1--10, then 11--20, then 21--30, ... if BATCH_INCREMENT = 10.

  // Declare some other preferences to be set later...
  long MULTIPLICTY, BATCH_MIN, BATCH_MAX;

  
  // Declare some empty vectors
  vector<long> empty_vec;
  vector<Matrix_mpz> empty_matrix_vec;
  vector<LocalConditions> empty_local_conditions_vec;


  /*
  // DIAGNOSTIC:
  cout << " Entering LocalSplitCoverings::_IncrementMinimalLocalSplitCoverings()" << endl; 
  cout << "   __Largest_size_used = " << __Largest_size_used << endl;
  cout << "   __used_indices = " << __used_indices << endl;
  cout << "   __complementary_matrix_vec = " << __complementary_matrix_vec << endl;
  cout << "   __local_conditions_vec = " << __local_conditions_vec << endl;
  */


  // Check if there is any previous temporary data
  if ( (__Largest_size_used == 0) || 
       (__complementary_matrix_vec.empty() == true) || (__local_conditions_vec.empty() == true) ) {

    // Make sure all temporary data is empty...
    assert(__Largest_size_used == 0);
    assert(__complementary_matrix_vec.empty() == true);
    assert(__local_conditions_vec.empty() == true);

    // Initialize the lists of complimentary matrices and their local conditions (to account for the zero vector).
    __complementary_matrix_vec.resize(1, empty_matrix_vec);
    __local_conditions_vec.resize(1, empty_local_conditions_vec);  


    // Set the desired multiplicity of the covering
    MULTIPLICTY = 1;
    
    // Set the size for making each batch of short vectors
    BATCH_MIN = 1;
    BATCH_MAX = BATCH_INCREMENT;
  }

  else {    

    // Make sure there are no (egregious) discrepancies in the temporary data
    assert(__Largest_size_used > 0);
    assert(__used_indices.empty() == false);
    /*
    cout << " __Largest_size_used = " << __Largest_size_used << endl; 
    cout << " __complementary_matrix_vec.size() = " << __complementary_matrix_vec.size() << endl; 
    cout << " __complementary_matrix_vec = " << __complementary_matrix_vec << endl; 
    */
    assert(__complementary_matrix_vec.size() == __Largest_size_used + 1);
    assert(__complementary_matrix_vec.size() == __local_conditions_vec.size());


    // Set the desired (additional) multiplicity to add one local cover
    MULTIPLICTY = 1;
    
    // Set the size for making each batch of short vectors
    BATCH_MIN = __Largest_size_used + 1;
    BATCH_MAX = __Largest_size_used + BATCH_INCREMENT;

  }




  // ---------------------- Now make the local covering by finding its indices --------------------------


  // Some local variables
  bool covering_flag = false;
  vector< set< vector<long> > > covering_indices;     // A list of sets of indices of the form (d,i)


  // Loop through batches of short vectors to look for local coverings
  while (covering_flag == false) {

    // Make the vectors with Q(v) <= BATCH_MAX
    vector< vector< vector<long> > > tmp_theta_vec;
    tmp_theta_vec = _ThetaVectors(_big_form, BATCH_MAX); 
    
    // Make an associated list of complementary matrices
    __complementary_matrix_vec.insert(__complementary_matrix_vec.end(), BATCH_INCREMENT, empty_matrix_vec); 
    for(long i=BATCH_MIN; i <= BATCH_MAX; i++)       // Note: The zero vector is ignored, since it doesn't give a complementary sublattice!
      for(long j=0; j < tmp_theta_vec[i].size(); j++)
	__complementary_matrix_vec[i].push_back(_FindComplementarySublattice(_big_form, tmp_theta_vec[i][j]));
    
    
    // Make an associated list of local conditions
    __local_conditions_vec.insert(__local_conditions_vec.end(), BATCH_INCREMENT, empty_local_conditions_vec);
    for(long i=BATCH_MIN; i <= BATCH_MAX; i++)       // Note: The zero vector is ignored, since it doesn't give a complementary sublattice!
      for(long j=0; j < tmp_theta_vec[i].size(); j++) {
	
	// CHOICE #1: Split quaternary local cover
	// ----------
	
	// This is for a split quaternary covering (EASIER TO FIND)
	Matrix_mpz split_matrix(4,4);
	split_matrix(1,1) = 2*i;  // Set the split value
	for(long ii=1; ii<=3; ii++)
	  for(long jj=1; jj<=3; jj++)
	    split_matrix(ii+1,jj+1) = __complementary_matrix_vec[i][j](ii,jj);
	
	// /*
	// DIAGNOSTIC
	cout << " About to find the local conditions for i = " << i << "  j = " << j << endl;    //" with the matrix vector: " << matrix_vec[i] << endl;
	//cout << " About to find the local conditions for " << endl << split_matrix << endl << endl;
	// */
	
	// Add the new local condition
	LocalConditions new_condition(split_matrix);   
	__local_conditions_vec[i].push_back(new_condition);      
      }
    
    
    // Try to construct MULTIPLICITY split local coverings
    covering_indices = _Find_Single_Form_Local_Covers(__local_conditions_vec, covering_indices, __used_indices);
    
    
    // /*
    // DIAGNOSTIC
    cout << endl << endl;
    cout << " The (cumulative) covering indices are " << covering_indices << endl;
    cout << " The _big_form local conditions are: " << _big_local << endl;
    cout << "==============================================================" << endl;
    cout << endl << endl;
    // */
    
    
    // Check if we're done!
    if (covering_indices.size() >= MULTIPLICTY)                   // Check if we have a covering
      covering_flag = true;
    else {
      BATCH_MIN = BATCH_MAX + 1;
      BATCH_MAX += BATCH_INCREMENT;
      //      exit(1);                                             ///////  PREMATURE EXIT....  DEAL WITH THE LOOP LATER...
    }

  }


  // -----------------------------  Contruct the forms from their indices ------------------------------------------

  // Sort out the original forms used from the covering indices 
  for(long ind=0; ind < covering_indices.size(); ind++) {

    // Find the new index we're adding
    long new_ind = _d_list.size();

    // Append an extra empty vector for each local cover
    _d_list.push_back(empty_vec);
    _form_list.push_back(empty_matrix_vec);
    _local_list.push_back(empty_local_conditions_vec);

    for(set< vector<long> >::iterator i = covering_indices[ind].begin(); i != covering_indices[ind].end(); i++) {
      
      // Add the correct information to the return lists
      _d_list[new_ind].push_back((*i)[0]);
      _form_list[new_ind].push_back(__complementary_matrix_vec[(*i)[0]][(*i)[1]]);
      _local_list[new_ind].push_back(__local_conditions_vec[(*i)[0]][(*i)[1]]);
      
    }
  }
  
  /*
  // Print all of the info
  cout << "Looking at the form: " << endl;
  cout << _big_local << endl;
  cout << " The index list is: " << endl;
  cout << index_list << endl;
  cout << " The form list is: " << endl;
  cout << form_list << endl;
  cout << " The d_list is: " << endl;
  cout << d_list << endl;
  cout << " The local conditions list is: " << endl;
  cout << local_list << endl;
  cout << " The running conditions list is: " << endl;
  cout << running_local_list << endl;
  cout << endl << endl << endl;
  */


  // Write the temporary data in case we need to add more covers later
  __Largest_size_used = BATCH_MAX;        // Says the length of the largest vector previously used.

  // Write the (now extended) local cover to a file (if we set an output filename)
  _WriteLocalSplitCoverings();
        
}






// ----------------------------------------------------------------------------------------------------------



// Checks if a vector of local coverings admits single form local covers
// (Returns empty set <==> No form is a local cover!)
vector< set< vector<long> > > LocalSplitCoverings::_Find_Single_Form_Local_Covers( 
										 const vector< vector<LocalConditions> > & local_vec, 
										 const vector< set< vector<long> > > & old_cover_indices,
										 set< vector<long> > & excluded_indices) const {
  

  
  // Make the vector of indexing sets to return
  vector< set< vector<long> > > return_list;
  return_list = old_cover_indices;


  // Some settings here
  long MULTIPLICITY = 5;                   // **********  THIS SHOULD BE PASSED IN !!!  ***********



  // 1. Look for local covers using only one form:
  // ---------------------------------------------

  // Run through all indices not in the excluded set
  vector<long> tmp_index;
  tmp_index.resize(2, 0);
  for(long i=0; i < local_vec.size(); i++)
    for(long j=0; j < local_vec[i].size(); j++) {
      tmp_index[0] = i;
      tmp_index[1] = j;

      // Check that the index isn't excluded, and we're not done.
      if ((return_list.size() < MULTIPLICITY) && (excluded_indices.find(tmp_index) == excluded_indices.end())) { 

	// Check if the form is a local cover by itself
	if (local_vec[i][j] == _big_local) {
	  set< vector<long> > univ_set;
	  univ_set.insert(tmp_index);	
	  return_list.push_back(univ_set);
	  excluded_indices.insert(tmp_index);
	}	
	
      }
    }
  
  // Check if we're done
  return return_list;

}



// ----------------------------------------------------------------------------------------------------------

/////////////////////////////////////////////////////
// Writes a local split covering to a datafile!
////////////////////////////////////////////////////
void LocalSplitCoverings::_WriteLocalSplitCoverings() const { 

  // Write the covering if we have a non-empty filename
  if(_filename_string.empty() == false) {
    
    // Try opening data file
    ofstream outfile;
    outfile.open(_filename_string.c_str());
    if (! outfile.is_open())
      { cout << "_WriteLocalSplitCoverings Error: Error opening file"; exit (1); }
    
    // Write and close the file
    _FileOutput(outfile); 
    outfile.close();
    
  }

}


/////////////////////////////////////////////////////////////////////////////////////////
// Reads a local split covering form a datafile, and checks it covers the original form!
////////////////////////////////////////////////////////////////////////////////////////
void LocalSplitCoverings::_ReadLocalSplitCoverings() {

  // Attempt to read the covering if we have a non-empty filename
  if(_filename_string.empty() == false) {
    
    // Check that the file exists
    if (FileExists(_filename_string) == false) {
      cout << "Error in _ReadLocalSplitCoverings: The filename " << _filename_string << " doesn't exist! =( " << endl;  
      assert( 1 == 0 );
    }
    
    // Try opening data file
    ifstream infile;
    infile.open(_filename_string.c_str());
    if (! infile.is_open())
      { cout << "_ReadLocalSplitCoverings Error: Error opening file"; exit (1); }
    
    // Read and close the file
    _FileInput(infile);   // Note: This includes some sanity checks! =)
    infile.close();
    
  }

}








// ============================================================================================================




// Test the ThetaVectorsRoutine by checking it against QQ.compute() by
// evaluating all of the vectors, checking they have the correct value
// and there are the right number of them.

void LocalSplitCoverings::_Testing_ThetaVectors() const {

  Matrix_mpz QQ;
  QQ.SetDims(4,4);
  /*
  // Using x^2 + 3y^2 + 5z^2 + 7w^2  (global matrix is 2Q)
  QQ(1,1) = 2;
  QQ(2,2) = 6;
  QQ(3,3) = 10;
  QQ(4,4) = 14;
  */

  // Using x^2 + 2y^2 + yz - yw + 4z^2 + 3zw + w^2  (global matrix is 2Q)
  QQ(1,1) = 2;
  QQ(2,2) = 4;
  QQ(2,3) = 1;
  QQ(3,2) = 1;
  QQ(2,4) = -1;
  QQ(4,2) = -1;
  QQ(3,3) = 8;
  QQ(3,4) = 3;
  QQ(4,3) = 3;
  QQ(4,4) = 62;



  unsigned long Precision = 100;  // This is the biggest we want the vectors

  cout << "Starting Testing_ThetaVectors()" << endl << endl;

  // Make the theta vectors
  cout << "Starting to make the theta vector" << endl;
  vector< vector< vector <long> > > Theta_vec;
  Theta_vec = _ThetaVectors(QQ, Precision);
  cout << "Finished making the theta vector" << endl;

  // Check the vectors have the correct length
  for(long i=1; i <= Precision; i++) {
    cout << " Checking i = " << i << endl;
    for(long j=0; j < Theta_vec[i].size(); j++){
      mpz_class value;
      value = QQ.EvaluateQuadratic(Theta_vec[i][j]);
      if (value != 2*i)                                               //  We multiply by 2 since we're using the global matrix of 2Q
	cout << "Error in Testing_ThetaVectors(): The value " << (value/2) 
	     << " of the vector " << Theta_vec[i][j] 
	     << " doesn't match its location " << i << endl;           // Note: We use i+1 since we don't count the zero vector!
      /*
      else
	cout << " Vector " << Theta_vec[i][j] << " is ok. =) " << endl;  
      */

    }
  }


  // Check the number of vectors agrees with the theta coefficient -- TO DO
  /*
  vector<long> Theta;
  QQ.compute()
  */


  cout << "Finishing Testing_ThetaVectors()" << endl;


}




// ================================ Strict I/O Front-end Routines ==========================================


///////////////////////////////////////////////
// Writes the local conditions to the ostream
///////////////////////////////////////////////
void FileOutput(const LocalSplitCoverings & elt, ostream & out) {

  elt._FileOutput(out);

}


////////////////////////////////////////////////
// Reads the local conditions from the istream
////////////////////////////////////////////////
void FileInput(LocalSplitCoverings & elt, istream & in) {

  elt._FileInput(in);

}


