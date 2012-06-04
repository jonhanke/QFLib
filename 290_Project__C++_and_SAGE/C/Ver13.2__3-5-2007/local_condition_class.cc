


// Empty Constructor
LocalConditions::LocalConditions(){
  _empty_flag = true;
  // Do nothing...
}

// Constructor using "Make(QQ)" to initialize the local conditions
LocalConditions::LocalConditions(const Matrix_mpz & QQ, const string & type){
  
  // Sanity Check
  assert(QQ.IsQuadraticForm() == true);
  assert((type == "slow_at_2") || (type == "fast_at_2"));


  // Choose the method
  if  (type == "slow_at_2")
    Make(QQ);
  else if (type == "fast_at_2")
    Make2(QQ);             //  <---- This one separates out the local computation in Make_Local_Repn_Vector_At_Prime()

}


// Empty Destructor 
LocalConditions::~LocalConditions(){
  // Do nothing...
}




vector<long> LocalConditions::Get_Anisotropic_Primes() const {
  if(_aniso__mod_meaning_flag == true)
    return aniso_vector;  
  else {
    cout << "Error in Get_Anisotropic_Primes():  The anisotropic vector is meaningless! =(" << endl;
    assert( 0 == 1 );
  }
};


vector<long> LocalConditions::Get_Local_Mod_Vector() const {
  if(_aniso__mod_meaning_flag == true)
    return local_mod_vector;  
  else {
    cout << "Error in Get_Local_Mod_Vector():  The anisotropic vector is meaningless! =(" << endl;
    assert( 0 == 1 );
  }
};


vector< vector<long> > LocalConditions::Get_Local_Repn_Array() const {
  return local_repn_array;
};




// Lists the primes with a congruence obstruction
vector<long> LocalConditions::List_Bad_Primes() const {

  vector<long> bad_primes;
  
  for(long i=0; i<local_repn_array.size(); i++) {
    
    bool has_obstruction = false;
    
    // Check each prime to see if it has an obstruction
    if (local_repn_array[i][0] == 2) {    // Check if p = 2 has an obstruction
      for(long j=1; j<9; j++)
	if (local_repn_array[i][j] != 0)
	  has_obstruction = true;
    } 
    else {                                 // Check if p > 2 has an obstruction
      for(long j=1; j<4; j++)
	if (local_repn_array[i][j] != 0)
	  has_obstruction = true;
    } 
    
    // Add obstructions to the list
    if (has_obstruction == true)
      bad_primes.push_back(local_repn_array[i][0]);

  }

  return bad_primes;

};






void LocalConditions::IntersectWith(const LocalConditions & B) {

  // Note: From this operation, we see that the anisotropic primes have no meaning here, 
  // and also the local mods...

  // Return the empty local condition if either of them is empty
  if (((*this).IsEmpty() == true) || (B.IsEmpty() == true))
    (*this)._empty_flag = true;
  
  else {
    
    // Make a vector of the union of these vectors
    vector< vector<long> > new_array;
    new_array = local_repn_array;
    for(long i=0; i<B.local_repn_array.size(); i++)
      new_array.push_back(B.local_repn_array[i]);
    
    
    // Sort by the primes to put them in order (Bubble sort)
    // Note: This could be more efficient since the vectors are increasing 
    // (i.e., the primes are in increasing order in each vector...)
    bool done_flag = false;
    vector<long> swap_vec;
    while (done_flag == false) {
      
      done_flag = true;
      
      for (long i=0; i+1 < new_array.size(); i++)
	if (new_array[i+1][0] < new_array[i][0]) {
	  done_flag = false;
	  
	  for(long j=0; j<9; j++)
	    swap_vec[j] = new_array[i][j];
	  for(long j=0; j<9; j++)
	    new_array[i][j] = new_array[i+1][j];
	  for(long j=0; j<9; j++)
	    new_array[i+1][j] = swap_vec[j];
	}
    }
    
    
    // If two adjacent ones are the same, then conglomerate them...
    for (long i=0; i+1 < new_array.size(); i++) 
      if (new_array[i+1][0] == new_array[i][0]) {
	
	// Make the new vector
	for(long j=0; j<9; j++) 
	  if ((new_array[i][j] >=0) && (new_array[i+1][j] >=0))
	    new_array[i][j] = max(new_array[i][j], new_array[i+1][j]);
	  else 
	    new_array[i][j] = -99;
	
	// Remove the duplicate entry
	new_array.erase(new_array.begin() + i+1);
	
      }
    
    
    // Replace the current entries with the intersection
    local_repn_array = new_array;
    aniso_vector.clear();
    local_mod_vector.clear();
    _aniso__mod_meaning_flag = false;
    
  }
  
}



void LocalConditions::UnionWith(const LocalConditions & B) {

  // Forget about this routine!
  cout << "ERROR: LocalConditions::UnionWith() should not be called! =0 " << endl;
  assert(0==1);

}




// This is used only in the
// LocalSplitCovering::_Find_Eligible_Covering_Primes(), since it does
// exactly what we want there for our necessary condition.
void LocalConditions::PrimewiseUnionWith(const LocalConditions & B) {

  // ==================================================================================

  // Note: From this operation, we see that the (private variables for
  // the) anisotropic primes have no meaning here, and also the local
  // mods...


  // Note: In contrast to the IntersectionWith() routine, here we must
  // check that each prime appears in at least once in each set of
  // local conditions to keep it!  To do this we assume that:
  //       * The primes are not repeated in any local condition
  //         (so 2 primes appearing is necessary to keep the condition!)
  //
  // To Do: For a speed improvement, we could also assume the primes
  // are in increasing order, which is probably already true, but I
  // haven't checked.




  // /*
  // DIAGNOSTIC
  cout << " \n\n Entering PrimewiseUnionWith(): " << endl;
  cout << " The number of new local conditions in B is: " << B.local_repn_array.size() << endl;
  cout << " Finding the union of the local conditions: \n " << B << endl;
  cout << " with the local conditions: \n " << (*this) << endl;
  // */





  // If both are empty, then the union is empty
  if (((*this).IsEmpty() == true) && (B.IsEmpty() == true))
    (*this)._empty_flag = true;
  
  
  
  // If (*this) is empty, then copy B into it
  else if ((*this).IsEmpty() == true) {
    (*this).local_repn_array = B.local_repn_array;
    (*this)._empty_flag = B._empty_flag;
  }


  // Otherwise, be sure that B is non-empty before continuing (so we have repeated conditions and don't kill everything unless B is full...)
  else if (B.IsEmpty() == false) {

    
    
    // Make a vector of the union of these vectors
    vector< vector<long> > new_array;
    new_array = local_repn_array;
    for(long i=0; i<B.local_repn_array.size(); i++)
      new_array.push_back(B.local_repn_array[i]);
    
    
    // /*
    // DIAGNOSTIC
    cout << " new_array.size() = " << new_array.size() << endl;
    cout << " new_array = " << endl << new_array << endl;
    // */
    
    
    
    // Sort by the primes to put them in order (Bubble sort)
    // Note: This could be more efficient since the vectors are increasing 
    // (i.e., the primes are in increasing order in each vector...)
    bool done_flag = false;
    vector<long> swap_vec(9,0);   // Starts out with length 9 and full of 0's
    while (done_flag == false) {
      
      done_flag = true;
      
      for (long i=0; i+1 < new_array.size(); i++)
	if (new_array[i+1][0] < new_array[i][0]) {
	  done_flag = false;
	  
	  for(long j=0; j<9; j++)
	    swap_vec[j] = new_array[i][j];
	  for(long j=0; j<9; j++)
	    new_array[i][j] = new_array[i+1][j];
	  for(long j=0; j<9; j++)
	    new_array[i+1][j] = swap_vec[j];
	}
    }
    
    
    // Remove Singleton Local Conditons:
    // ---------------------------------
    
    // Make sure only to keep repeated primes.         // Note: We assume that the primes only appear with multiplicity 1 or 2!
    unsigned long ptr = 0;
    while (ptr + 1 < new_array.size())
      if (new_array[ptr+1][0] != new_array[ptr][0]) {
	new_array.erase(new_array.begin() + ptr);
      }
      else
	ptr = ptr + 2;
    
    // If we're at then end, then that local condition appears by itself (by the parity of our construction)
    if (ptr + 1 == new_array.size())
      new_array.erase(new_array.begin() + ptr);
    
    
    
    
    
    // /*
    // DIAGNOSTIC
    cout << " new_array.size() = " << new_array.size() << endl;
    cout << " new_array.size() - 1 = " << (new_array.size() - 1) << endl;
    cout << " new_array = " << endl << new_array << endl;
    // */
    
    
    
    // Conglomerate Duplicate Local Conditons:
    // ---------------------------------------
    
    // If two adjacent ones are the same, then conglomerate them...
    for (long i=0; i+1 < new_array.size(); i++) 
      if (new_array[i+1][0] == new_array[i][0]) {
	
	// Make the new vector
	for(long j=0; j<9; j++) 
	  if ((new_array[i][j] >=0) && (new_array[i+1][j] >=0))         // If both are >=0 take the smallest
	    new_array[i][j] = min(new_array[i][j], new_array[i+1][j]);
	  else if ((new_array[i][j] >=0) || (new_array[i+1][j] >=0))    // If only one is >=0, then take the biggest
	    new_array[i][j] = max(new_array[i][j], new_array[i+1][j]);
	  else                                                          // Otherwise, set it to -99
	    new_array[i][j] = -99;
	
	// Remove the duplicate entry
	new_array.erase(new_array.begin() + i+1);
	
      }
    
    
    // Replace the current entries with the union
    local_repn_array = new_array;
    aniso_vector.clear();
    local_mod_vector.clear();
    _aniso__mod_meaning_flag = false;
    
    // Make sure that the existing entries are all meaningful
    // (i.e., it could happen that the union of two partial obstructons gives no obstruction.)
    Cleanup();
    
  }


  // /*
  // DIAGNOSTIC
  cout << " Exiting UnionWith() with the local conditions: \n " << (*this) << endl << endl;
  // */
}












// Get the local conditions for a (global) quadratic form 
// -- This is just the "FindLocalConditions" routine from new_local_representability.cc
void LocalConditions::Make(const Matrix_mpz & QQ) {

  // Wipe any old information here!
  aniso_vector.clear();
  local_mod_vector.clear();
  local_repn_array.clear();
  _empty_flag = true;
  _aniso__mod_meaning_flag = true;

  // WARNING:  This only works when the number of level primes is <= the number of rows in local_repn_array...   

/*
void FindTernaryLocalConditions(long local_array[][9], long Q[6],   
vector<long> & prime_list) { 
*/

  // NOTATION: 
  // =========
  // The structure of the local array is given by:
  //         [ p, *, *, *, *, 0, 0, 0, 0]  for p>2
  //   where the * is a number saying the smallest power k so that the
  //   number (p^2k * x) in the squareclass x is represented.  When
  //   p>2 this is given by x = 1, r, p, p*r.  When p=2 this is given
  //   by x = ..... .


  // Functions needed:
  //   ord(x, p)
  //   QFlevel(Q) -- long QF_Ternary_Level(long[6])
  

  /*
  // DIAGNOSTIC
  // Show the matrix QQ
  // (By convention QQ will be the matrix of 2Q for integrality reasons)
  cout << " QQ is given by: " << endl;
  cout << QQ << endl;
  */


  // Set the Local Condtions to be non-empty (since no form has empty local conditions...
  (*this)._empty_flag = false;


  // Find the level
  //      long N = QF_Ternary_Level(Q); 
  mpz_class NN = QQ.QFLevel();
  //  long N;  N = NN.get_ui();

  // DIAGNOSTIC



  /*
  // Remove any overall scaling factor (to ensure QQ is primitive)
  mpz_class scaling_factor = 1;
  for(i=1; i<QQ.NumRows(); i++)
    for(j=i; j<QQ.NumRows(); j++)
      if (i==j)
	scaling_factor = GCD(scaling_factor, QQ(i,i));
      else
	scaling_factor = GCD(2*scaling_factor, QQ(i,i));
  */




  // Find the number of primes dividing the level
  vector<long> level_primes = PrimeDivisors1(2*NN);
  long Num_Primes = level_primes.size();


  /*
  // DIAGNOSTIC
  cout << " The level is NN = " << NN << endl;
  cout << " Num_Primes: " << Num_Primes << endl;
  cout << " Level Primes: " << level_primes << endl;
  */
  

  // Find the modulus for each prime 
  long repn_mods[Num_Primes];
  for (long i=0; i<Num_Primes; i++)
    repn_mods[i] = LongPow(level_primes[i], Valuation(4*NN, level_primes[i]) + 2);
 

  /*
  // DIAGNOSTIC: PRINT THE LOCAL CONDITION INFO...
  cout << "NN = " << NN << endl;
  for (long i=0; i<Num_Primes; i++)
    cout << " i = " << i << "   level prime = " << level_primes[i] 
	 << "   repn_mod = " << repn_mods[i] << endl;
  */


  // Make a table of local normal forms of QQ at each prime p | N
  Matrix_mpz local_normal_forms[Num_Primes];
  for (long i=0; i < Num_Primes; i++)
    local_normal_forms[i] = QQ.GetLocalNormal(mpz_class(level_primes[i]));

  /*
  // Make the representation array
    long repn_array[Num_Primes][9];
  */


  /*
  // DIAGNOSTIC
  cout << " Here is the uninitialized array " << endl;
  cout << " These are the local representability conditions (array): " << endl;
  cout << local_repn_array << endl;
  cout << endl;
  */


  // Check local representability for each prime
  local_repn_array.clear();

  for (long i=0; i<Num_Primes; i++) {
    vector<long> tmp_local_repn_vec(9, 0);   // Temporary vector for the local conditions at a fixed prime, initialized to all zeros!    
    mpz_class p = level_primes[i];
    tmp_local_repn_vec[0] = p.get_ui();
    
    // Make the squareclass vector
    long sqclass_size;
    if (p == 2) 
      sqclass_size = 8;
    else 
      sqclass_size = 4;

    long sqclass[sqclass_size];

    if (p == 2) {
      sqclass[0] = 1;
      sqclass[1] = 3;
      sqclass[2] = 5;
      sqclass[3] = 7;
      sqclass[4] = 2;
      sqclass[5] = 6;
      sqclass[6] = 10;
      sqclass[7] = 14;
    }
    else {
      long r = NonResidue(p).get_ui();
      sqclass[0] = 1;
      sqclass[1] = r;
      sqclass[2] = p.get_ui();
      sqclass[3] = p.get_ui()*r;
    }

    // Check the representability in each squareclass
    for (long j=0; j<sqclass_size; j++) {
      mpz_class m = sqclass[j];    // We need to use an mpz_class here to ensure there's no overflow! =)
      long k = 0;
      bool repn_flag = false;
      //      cout << "hi -- Using p = " << level_primes[i] << " and 4*NN = " << (4*NN) << endl;
      long prime_pow = Valuation(4*NN, level_primes[i]) + 2;  // The power of the repn modulus
      //      cout << "bye" << endl;
      while ((repn_flag == false) && (m < 4*NN*p*p)) {    // Note: We used the mpz_class NN here just to be safe.

	/*
	// DIAGNOSTIC
	cout << endl;
	cout << " m = " << m << endl;
	cout << " p = " << p << endl;
	cout << " Local Normal Form = " << endl << local_normal_forms[i];
	cout << " The Local Density is " << local_normal_forms[i].Local_Density(p, m) << endl;
	*/

	if (local_normal_forms[i].Local_Density(p, m) > 0) {
	  tmp_local_repn_vec[j+1] = k;
	  repn_flag = true;
	}
	k++;
	m = m * p * p;
      }
      
      // If we're not represented, write a negative number
      // to signify we checked up to (but not including) x * p^(2*k).
      if (repn_flag == false)
	tmp_local_repn_vec[j+1] = -k;
      
    }

    local_repn_array.push_back(tmp_local_repn_vec);

  }

  /*
  // DIAGNOSTIC
  cout << " Here is the initialized array " << endl;
  cout << " These are the local representability conditions (array): " << endl;
  cout << local_repn_array << endl;
  cout << endl;
  */

      
  // Make the big modulus
  long Bigmod = 1;
  for (long i=0; i < Num_Primes; i++)
    Bigmod = Bigmod * repn_mods[i];


  // Return the modulus and list of local failures
  //       local_repn_array = repn_array;
  
  local_mod_vector.clear();
  for(long j=0; j<Num_Primes; j++)
    local_mod_vector.push_back(repn_mods[j]);

  /*
  // DIAGNOSTIC
  cout << "\n\n Making Aniso Primes \n\n";
  */


  aniso_vector.clear();
  for(long j=0; j<Num_Primes; j++) {

    /*
    // DIAGNOSTIC
    cout << " ... Using the form " << endl << QQ;
    cout << " IsAnisotropic(QQ,2) = " << QQ.IsAnisotropic(2) << endl;
    cout << " IsAnisotropic(QQ,29) = " << QQ.IsAnisotropic(29) << endl;
    */

    if (QQ.IsAnisotropic(mpz_class(level_primes[j])) == true) {
      aniso_vector.push_back(level_primes[j]);

      /*
      // DIAGNOSTIC
      cout << " Found an anisotropic prime: " << level_primes[j] << endl;
      */
    }
  }

  //   cout << "\n\n Finished Making Aniso Primes \n\n";


  // Clean up the local_representation conditions
  Cleanup();

}





// Get the local conditions for a (global) quadratic form 
// -- This is just the "FindLocalConditions" routine from new_local_representability.cc
void LocalConditions::Make2(const Matrix_mpz & QQ) {

  // Wipe any old information here!
  aniso_vector.clear();
  local_mod_vector.clear();
  local_repn_array.clear();
  _empty_flag = true;
  _aniso__mod_meaning_flag = true;

  // NOTATION: 
  // =========
  // The structure of the local array is given by:
  //         [ p, *, *, *, *, 0, 0, 0, 0]  for p>2
  //   where the * is a number saying the smallest power k so that the
  //   number (p^2k * x) in the squareclass x is represented.  When
  //   p>2 this is given by x = 1, r, p, p*r.  When p=2 this is given
  //   by x = ..... .

  

  /*
  // DIAGNOSTIC
  // Show the matrix QQ
  // (By convention QQ will be the matrix of 2Q for integrality reasons)
  cout << " QQ is given by: " << endl;
  cout << QQ << endl;
  */


  // Set the Local Condtions to be non-empty (since no form has empty local conditions...
  (*this)._empty_flag = false;


  // Find the level
  mpz_class NN = QQ.QFLevel();


  // Find the number of primes dividing the level
  vector<long> level_primes = PrimeDivisors1(2*NN);
  long Num_Primes = level_primes.size();


  /*
  // DIAGNOSTIC
  cout << " The level is NN = " << NN << endl;
  cout << " Num_Primes: " << Num_Primes << endl;
  cout << " Level Primes: " << level_primes << endl;
  */
  

  // Find the modulus for each prime 
  long repn_mods[Num_Primes];
  for (long i=0; i<Num_Primes; i++)
    repn_mods[i] = LongPow(level_primes[i], Valuation(4*NN, level_primes[i]) + 2);
 

  /*
  // DIAGNOSTIC: PRINT THE LOCAL CONDITION INFO...
  cout << "NN = " << NN << endl;
  for (long i=0; i<Num_Primes; i++)
    cout << " i = " << i << "   level prime = " << level_primes[i] 
	 << "   repn_mod = " << repn_mods[i] << endl;
  */


  // Make a table of local normal forms of QQ at each prime p | N
  Matrix_mpz local_normal_forms[Num_Primes];
  for (long i=0; i < Num_Primes; i++)
    local_normal_forms[i] = QQ.GetLocalNormal(mpz_class(level_primes[i]));


  // Check local representability for each prime
  local_repn_array.clear();
  for (long i=0; i<Num_Primes; i++) 
    local_repn_array.push_back(Make_Local_Repn_Vector_At_Prime(local_normal_forms[i], NN, level_primes[i]));

  

  /*
  // DIAGNOSTIC
  cout << " Here is the initialized array " << endl;
  cout << " These are the local representability conditions (array): " << endl;
  cout << local_repn_array << endl;
  cout << endl;
  */

      
  // Make the big modulus
  long Bigmod = 1;
  for (long i=0; i < Num_Primes; i++)
    Bigmod = Bigmod * repn_mods[i];


  // Return the modulus and list of local failures
  //       local_repn_array = repn_array;
  
  local_mod_vector.clear();
  for(long j=0; j<Num_Primes; j++)
    local_mod_vector.push_back(repn_mods[j]);

  /*
  // DIAGNOSTIC
  cout << "\n\n Making Aniso Primes \n\n";
  */

  
  // Find the anisotropic primes
  aniso_vector.clear();
  for(long j=0; j<Num_Primes; j++) {
    
    if (QQ.IsAnisotropic(mpz_class(level_primes[j])) == true) {
      aniso_vector.push_back(level_primes[j]);
      
      /*
      // DIAGNOSTIC
      cout << " Found an anisotropic prime: " << level_primes[j] << endl;
      */
    }
  }

  //   cout << "\n\n Finished Making Aniso Primes \n\n";


  // Clean up the local_representation conditions
  Cleanup();

}





// Get the local conditions for a local normalized quadratic form at a prime p 
vector<long> LocalConditions::Make_Local_Repn_Vector_At_Prime(const Matrix_mpz & local_normal_at_p, const mpz_class & NN, const mpz_class & p) const {

  // NOTATION: 
  // =========
  // The structure of the local array is given by:
  //         [ p, *, *, *, *, 0, 0, 0, 0]  for p>2
  //   where the * is a number saying the smallest power k so that the
  //   number (p^2k * x) in the squareclass x is represented.  When
  //   p>2 this is given by x = 1, r, p, p*r.  When p=2 this is given
  //   by x = ..... .
  

  // Sanity Checks:
  assert(local_normal_at_p.NumRows() >= 2);    // Check it's at least 2-dim'l.
  

  vector<long> local_conditions_at_p(9,0);
  local_conditions_at_p[0] = p.get_ui();



  // If p isn't a level prime, then it's universal at p! =)
  if ((NN % p != 0) && (p != 2)) {
    cout << "Warning: " << p << " is not a level prime..." << endl;
    return local_conditions_at_p;
  }


  
  //  mpz_class NN = QQ.QFLevel();
  
  
  
  // Find the modulus for each prime 
  mpz_class repn_mod = p ^ (Valuation(4*NN, p) + 2);
  
  
  /*
  // DIAGNOSTIC: PRINT THE LOCAL CONDITION INFO...
  cout << " Diagnostic for LocalConditions::Make_Local_Repn_Vector_At_Prime: " << endl;
  cout << " ---------------------------------------------------------------- " << endl;
  cout << " p = " << p << endl;
  cout << " Local Normal at p is " << endl << local_normal_at_p << endl;
  cout << " NN = " << NN << endl;
  cout << " repn_mod = " << repn_mod << endl;
  cout << endl;
  */
  
  
  
  
  // Make the squareclass vector
  long sqclass_size;
  if (p == 2) 
    sqclass_size = 8;
  else 
    sqclass_size = 4;
  
  long sqclass[sqclass_size];
  
  if (p == 2) {
    sqclass[0] = 1;
    sqclass[1] = 3;
    sqclass[2] = 5;
    sqclass[3] = 7;
    sqclass[4] = 2;
    sqclass[5] = 6;
    sqclass[6] = 10;
    sqclass[7] = 14;
  }
  else {
    long r = NonResidue(p).get_ui();
    sqclass[0] = 1;
    sqclass[1] = r;
    sqclass[2] = p.get_ui();
    sqclass[3] = p.get_ui()*r;
  }

  
  // Check the representability in each squareclass
  for (long j=0; j<sqclass_size; j++) {
    mpz_class m = sqclass[j];    // We need to use an mpz_class here to ensure there's no overflow! =)
    long k = 0;
    bool repn_flag = false;
    //      cout << "hi -- Using p = " << level_primes[i] << " and 4*NN = " << (4*NN) << endl;
    //    long prime_pow = Valuation(4*NN, p) + 2;  // The power of the repn modulus
    //      cout << "bye" << endl;
    while ((repn_flag == false) && (m < repn_mod)) {    // Note: We used the mpz_class NN here just to be safe.
      
      /*
      // DIAGNOSTIC
      cout << endl;
      cout << " m = " << m << endl;
      cout << " p = " << p << endl;
      cout << " Local Normal Form = " << endl << local_normal_at_p;
      cout << " The Local Density is " << local_normal_at_p.Local_Density(p, m) << endl;
      */
      
      if (local_normal_at_p.Local_Density(p, m) > 0) {
	local_conditions_at_p[j+1] = k;
	repn_flag = true;
      }
      k++;
      m = m * p * p;
    }
    
    // If we're not represented, write a negative number
    // to signify we checked up to (but not including) x * p^(2*k).
    if (repn_flag == false)
      local_conditions_at_p[j+1] = -k;
    
  }
  

  /*
  // DIAGNOSTIC
  cout << " The local conditions vector is " << local_conditions_at_p << endl;
  cout << endl;
  */


  // Return the local conditions at p
  return local_conditions_at_p;
  
}





void LocalConditions::Cleanup() {            // <---------  Question: Should we test this to see if it's empty?
  
  // Make a new local_repn_array
  vector< vector<long> > new_repn_array;

  // Add congruence obstruction "vectors" to it
  for(long i=0; i<local_repn_array.size(); i++) {

    bool has_obstruction = false;
    
    // Check each prime to see if it has an obstruction
    if (local_repn_array[i][0] == 2) {    // Check if p = 2 has an obstruction
      for(long j=1; j<=8; j++)
	if (local_repn_array[i][j] != 0)
	  has_obstruction = true;
    } 
    else {                                 // Check if p > 2 has an obstruction
      for(long j=1; j<=4; j++)
	if (local_repn_array[i][j] != 0)
	  has_obstruction = true;
    } 
    
    // Add obstructions to the list
    if (has_obstruction == true) {
      new_repn_array.push_back(local_repn_array[i]);
      //      cout << " Cleanup() Added : " << local_repn_array[i] << endl;
    }
  }

  // Now replace the local_repn_array with the cleaner one. =)
  local_repn_array = new_repn_array;

}





///////////////////////////////////////////////////////
// Compares when two local_repn_arrays are the same.
////////////////////////////////////////////////////////
bool LocalConditions::operator==(const LocalConditions & B) const {

  // Check they have the same size
  if (local_repn_array.size() != B.local_repn_array.size())
    return false;

  // Check they have the same entries
  for(long i=0; i<local_repn_array.size(); i++)
    if (local_repn_array[i] != B.local_repn_array[i]) 
      for (long j=0; j<=8; j++)
	if ((local_repn_array[i][j] != B.local_repn_array[i][j]) 
	    && ((local_repn_array[i][j] >= 0) || (B.local_repn_array[i][j] >= 0))) 
	  return false;

  // Check that they have the same empty flag
  if (_empty_flag != B._empty_flag)
      return false;


  // If both, then they're equal! =)
  return true;
}



///////////////////////////////////////////////////////
// Compares when two local_repn_arrays are different.
////////////////////////////////////////////////////////
bool LocalConditions::operator!=(const LocalConditions & B) const {

  return !((*this) == B);

}



///////////////////////////////////////////////////////
// Says if the local conditions are empty
////////////////////////////////////////////////////////
bool LocalConditions::IsEmpty() const {

  return _empty_flag;

}



///////////////////////////////////////////////////////
// Says if the local conditions are universal
////////////////////////////////////////////////////////
bool LocalConditions::IsUniversal() const {
  
  // Check it's not empty
  if (_empty_flag == true)
    return false;

  // Check it has no local obstructions
  if (local_repn_array.size() > 0)
    return false;

  
  return true;

}



///////////////////////////////////////////////////////////////
// Says if the local conditions are universal at the prime p
///////////////////////////////////////////////////////////////
bool LocalConditions::IsUniversalAtPrime(const mpz_class & p) const {

  // Sanity Check
  assert(p >= 2);

  
  // Check it's not empty
  if (_empty_flag == true)
    return false;

  // Check if it has no local obstructions
  if (local_repn_array.size() == 0)
    return true;

  // Check if the prime p appears
  for(long i=0; i<=local_repn_array.size(); i++) 
    if(local_repn_array[i][0] == p) {
      for(long j=1; j<=8; j++)
	if(local_repn_array[i][j] != 0)     // If we have some obstruction, return false!
	  return false;
    }
  
  // There are no obstructinos involving the prime p, so it's locally universal!
  return true;

}





/////////////////////////////////////////////////////
// Says whether or not a square-free m is locally represented.
// Requires m to be non-negative... for now.
/////////////////////////////////////////////////////

bool LocalConditions::IsLocallyRepresented__Squarefree(const mpz_class & m) const {

  // If the local conditions are empty, then m isn't represented
  if (_empty_flag == true)
    return false;

  // Check that m is positive.
  if (m < 0)
    return false;
  else if (m == 0)
    return true;


  // Otherwise, look for primes where it's not represented
  for(long i=0; i<local_repn_array.size(); i++) {

    long index;  // The index (between 1 and 8) of the local condition for m.
    long entry;  // The p-divisibility number associated to m in that index.
    
    // Find the entry corresponding to m.

    // Find the valuation and unit parts of m
    long p = local_repn_array[i][0];
    unsigned long v = Valuation(m, p);

    mpz_class m_unit;
    m_unit = m / (mpz_class(p)^v);     // Is this p^v notation ok?


    // Find the index of the unit part:
    // --------------------------------
    if (p == 2) {

      // If p = 2, reduce m_unit mod 8 and find the index
      m_unit = ((m_unit % 8) + 8) % 8;
      index = (m_unit.get_ui() + 1) / 2;

      // If it's an odd 2-power shift the index by 4
      if (v % 2 != 0)
	index += 4;
      
    } 
    else {
      
      // If p > 2, use the Legendre Symbol to find the index
      if (LegendreSymbol(m_unit, p) == 1)
	index = 1;
      else
	index = 2;
      
      // If it's an odd p-power shift the index by 2
      if (v % 2 != 0)
	index += 2;
      
    }
    
    
    // Find the entry for that index:
    // ------------------------------
    if (v % 2 == 0)      
      entry = v / 2;
    else
      entry = (v-1) / 2;


    // Return false if it's not locally represented at this prime:
    // -----------------------------------------------------------
    if ( (local_repn_array[i][index] < 0) || (entry < local_repn_array[i][index]) )
      return false;

  }



  // Since there are no local obstructions, it's represented! =)
  return true;

}



/////////////////////////////////////////////////////
// Says whether or not m is locally represented.
// Requires m to be non-negative... for now.
/////////////////////////////////////////////////////

bool LocalConditions::IsLocallyRepresented(const mpz_class & m) const {

  // If the local conditions are empty, then m isn't represented
  if (_empty_flag == true)
    return false;

  // Check that m is positive.
  if (m < 0)
    return false;
  else if (m == 0)
    return true;


  // Otherwise, look for primes where it's not represented
  for(long i=0; i<local_repn_array.size(); i++) {

    long index;  // The index (between 1 and 8) of the local condition for m.
    long entry;  // The p-divisibility number associated to m in that index.
    
    // Find the entry corresponding to m.

    // Find the valuation and unit parts of m
    long p = local_repn_array[i][0];
    unsigned long v = Valuation(m, p);

    mpz_class m_unit;
    m_unit = m / (mpz_class(p)^v);     // Is this p^v notation ok?


    // Find the index of the unit part:
    // --------------------------------
    if (p == 2) {

      // If p = 2, reduce m_unit mod 8 and find the index
      m_unit = ((m_unit % 8) + 8) % 8;
      index = (m_unit.get_ui() + 1) / 2;

      // If it's an odd 2-power shift the index by 4
      if (v % 2 != 0)
	index += 4;
      
    } 
    else {
      
      // If p > 2, use the Legendre Symbol to find the index
      if (LegendreSymbol(m_unit, p) == 1)
	index = 1;
      else
	index = 2;
      
      // If it's an odd p-power shift the index by 2
      if (v % 2 != 0)
	index += 2;
      
    }
    
    
    // Find the entry for that index:
    // ------------------------------
    if (v % 2 == 0)      
      entry = v / 2;
    else
      entry = (v-1) / 2;


    // Return false if it's not locally represented at this prime:
    // -----------------------------------------------------------
    if ( (local_repn_array[i][index] < 0) || (entry < local_repn_array[i][index]) )
      return false;

  }



  // Since there are no local obstructions, it's represented! =)
  return true;

}






////////////////////////////////
// Prints the local conditions
///////////////////////////////
ostream & LocalConditions::Print(ostream & out) const {    

  if ((*this).IsEmpty() == true)
    out << "Empty";
  else
    out << local_repn_array;

  return out;
}


////////////////////////////////////////////////////
// This is the overloaded (global) version for <<
////////////////////////////////////////////////////
ostream & operator<<(ostream & out, const LocalConditions & local) {    
  return local.Print(out);
}








///////////////////////////////////////////////
// Writes the local conditions to the ostream
///////////////////////////////////////////////
void LocalConditions::_FileOutput(ostream & out) const {

  // Opening Brace
  out << " << ";

  // Write the empty flag
  out << _empty_flag << " ; ";

  // Write the local modulus vector
  FileOutput(local_mod_vector, out);
  out << " ; ";

  // Write the anisotropic vector
  FileOutput(aniso_vector, out);
  out << " ; " << endl;

  // Write the local representation array
  FileOutput(local_repn_array, out);
  out << endl;

  // Closing Brace
  out << " >> " << endl;

}



////////////////////////////////////////////////
// Reads the local conditions from the istream
////////////////////////////////////////////////
void LocalConditions::_FileInput(istream & in) {

  // Declare temporary variables
  bool tmp_flag;
  vector<long> tmp_aniso_vec, tmp_local_mod_vec;
  vector< vector<long> > tmp_repn_array;
  char ch;
  long num;


  // Read the opening " << "
  in >> ch;
  assert(ch == '<');
  in >> ch;
  assert(ch == '<');

  // Read the empty_flag and ';'
  in >> tmp_flag;
  in >> ch;
  assert(ch == ';');

  // Read the local modulus vector and ';'
  FileInput(tmp_local_mod_vec, in);
  in >> ch;
  assert(ch == ';');

  // Read the anisotropic vector and ';'
  FileInput(tmp_aniso_vec, in);
  in >> ch;
  assert(ch == ';');

  // Read the local representation array
  FileInput(tmp_repn_array, in);

  // Read the closing " >> "
  in >> ch;
  assert(ch == '>');
  in >> ch;
  assert(ch == '>');


  // Perform some sanity checks 
  // TO DO: Check that all primes divide the level (which we must compute!) 


  // Copy the data to (*this)
  _empty_flag = tmp_flag;
  local_mod_vector = tmp_local_mod_vec;
  aniso_vector = tmp_aniso_vec;
  local_repn_array = tmp_repn_array;

}


// -----------------------------------------------------------------------------------------------

//////////////////////////////////////////////////////////
// Multiply a vector of local conditions at p=2 by 2^nu
//////////////////////////////////////////////////////////
vector<long> LocalConditions::_Multiply_Local_Conditions_by_power_at_p(const vector<long> & conditions, const long & nu, const mpz_class & p) const {

  // Sanity Checks:
  // --------------
  assert(p >= 2);
  assert(conditions.size() == 9);         // Check it represents a local conditions vector,
  assert(conditions[0] == p.get_ui());    // for the prime p = 2.
  assert(nu >= 0);                        // Check that the power 2^nu makes sense.


  // Declare the new conditions
  vector<long> new_conditions(9,0);
  new_conditions[0] = conditions[0];


  // Determine the shift
  long shift = nu/2;     // Shift each condition up by floor(nu/2);
  assert(shift + (nu % 2) == nu);

   

  // If p is 2...
  if (p==2) {
    
    // Deal with the parity of nu
    if (nu % 2 == 0)
      
      // Do nothing if nu is even
      new_conditions = conditions;    
    
    else {
      // If nu is odd
      for(long i=1; i<=4; i++) {
	
	// Reverse the squareclasses
	new_conditions[i] = conditions[i+4];      
	new_conditions[i+4] = conditions[i];            
	
	// Shift the (new) unit ones by 1, if it was represented
	if (new_conditions[i] >= 0)
	  new_conditions[i]++;
      }
    }
    
    // Finally, shift all the conditions up
    for(long i=1; i<=8; i++)
      new_conditions[i] += shift;
    
  }

  // Now if p is odd...
  else {
    
    
    // Determine the shift
    long shift = nu/2;     // Shift each condition up by floor(nu/2);
    assert(shift + (nu % 2) == nu);
    
    
    // Deal with the parity of nu
    if (nu % 2 == 0)
      
      // Do nothing if nu is even
      new_conditions = conditions;    
    
    else {
      // If nu is odd
      for(long i=1; i<=2; i++) {
	
	// Reverse the squareclasses
	new_conditions[i] = conditions[i+2];      
	new_conditions[i+2] = conditions[i];            
	
	// Shift the (new) unit ones by 1, if it was represented
	if (new_conditions[i] >= 0)
	  new_conditions[i]++;
      }
    }
    
    
    // Finally, shift all the conditions up
    for(long i=1; i<=4; i++)
      new_conditions[i] += shift;    
    
  }
  

  
  // Return the new conditions
  return new_conditions;
  
}




///////////////////////////////////////////////////////////////
// Multiply a vector of local conditions at p=2 by a unit u.
///////////////////////////////////////////////////////////////
vector<long> LocalConditions::_Multiply_Local_Conditions_by_unit_at_p(const vector<long> & conditions, const mpz_class & u, const mpz_class & p) const {

  // Sanity Checks:
  // --------------
  assert(p >= 2);
  assert(conditions.size() == 9);        // Check it represents a local conditions vector,
  assert(conditions[0] == p.get_ui());   // for the prime p = 2.
  assert(u % p != 0);                    // Check that u is a unit!


  if (p==2) {

    // Declare the new_conditions to return
    vector<long> new_conditions;
    new_conditions = conditions;    
    
    // Normalize the unit to its smallest positive squareclass representative
    long uu = (((u.get_ui() % 8) + 8) % 8);
    
    
    // Arrange for uu = 1 or 3 (mod 8)
    if ((uu == 5) || (uu == 7)) {
      
      // Reverse the squareclasses (u <=> -u) and continue
      uu = 8 - uu;
      for(long i=1; i<=2; i++) {      
	long tmp1 = new_conditions[i];
	new_conditions[i] = new_conditions[5-i];      
	new_conditions[5-i] = tmp1;
	
	long tmp2 = new_conditions[i+4];
	new_conditions[i+4] = new_conditions[9-i];            
	new_conditions[9-i] = tmp2;            
      }
      
    }
    
    
    // Finish the multplication by uu, and return the local conditions
    if (uu == 1)
      return new_conditions;
    else if (uu == 3) {
      for(long i=1; i<=4; i++) {      
	long tmp1 = new_conditions[2*i-1];
	new_conditions[2*i-1] = new_conditions[2*i];      
	new_conditions[2*i] = tmp1;
      }
      
    }

    else {
      
      if (LegendreSymbol(uu, p) == 1)
	return conditions;            // Do nothing if u is a square (mod p)
      else {

	// Declare the new_conditions to return
	vector<long> new_conditions(9,0);
	new_conditions[0] = conditions[0];
	new_conditions[1] = conditions[2];
	new_conditions[2] = conditions[1];
	new_conditions[3] = conditions[4];
	new_conditions[4] = conditions[3];
	return new_conditions;
      }
      
    }

  
    return new_conditions;
  }
  
}
    



////////////////////////////////////////////////////////////////////////////////////
// Parse a local normal form at some prime p into a vector of pairs describing it.
////////////////////////////////////////////////////////////////////////////////////
vector< vector<long> > LocalConditions::_Parse_Local_Normal_Form_at_p(const Matrix_mpz & normal_form, const mpz_class & p) const {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Code for pairs:
  // ---------------
  //   p > 2  ==>  (nu_i, 1 or -1)                          <-- Sign depends on the Legendre symbol
  //   p = 2  ==>  (nu_i, 1 or 3 or 5 or 7 or -1 or -3)
  //                  1,3,5,7 <--> x^2, 3x^2, 5x^2, 7x^2
  //                 -1 <--------> xy
  //                 -3 <--------> x^2 + xy + y^2
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  // Sanity Check
  assert(p > 2);

  // Test to see if it really is a local normal form at p.
  // ... DO THIS ...


  // Make the return vector
  vector< vector<long> > local_normal_vec;


  // Deal with prime p > 2 first
  if (p > 2) {
    mpz_class num;
    assert(num != 0);
    unsigned long pow;
    vector<long> pair(2,0);

    for(long i=1; i<=normal_form.NumRows(); i++) {
      // Read the 1-dim'l form 
      num = normal_form(i,i);
      pow = Valuation(num, p.get_ui());
      num = num / (p ^ pow);

      // Add the assocaited pair to the list
      pair[0] = pow;
      pair[1] = LegendreSymbol(num,p);
      local_normal_vec.push_back(pair);
    }

  }


  // Otherwise, p = 2
  else {
    mpz_class diag, off_diag;
    long pow;
    vector<long> pair(2,0);
    
    long i=1;
    while(i<=normal_form.NumRows()) {

      // Read the diagonal and off-diagonal coefficients
      diag = normal_form(i,i);
      if (i < normal_form.NumRows())
	off_diag = normal_form(i,i+1);     // Note: Here the matrix is a local notmal form, so it's the upper triangular matrix of Q!!!
      else 
	off_diag = 0;


      // ******** Sanity Check (to check that our normal form is triangular int he correct way...)
      assert((diag != 0) || (off_diag != 0));

	
      // Compute the 2-divisibility 
      mpz_class num;
      unsigned long pow;
      if (diag != 0)     
	num = diag;        // This deals with the forms a*x^2 and x^2 + xy + y^2.
      else 
	num = off_diag;    // This deals with the form x*y.
      pow = Valuation(num, p.get_ui());
      num = num / (p ^ pow);


      // Check if the form is 1-dim'l
      if (off_diag == 0) {
	
	// Add the associated pair to the list
	pair[0] = pow;
	pair[1] = (((num.get_ui() % 8) + 8) % 8);
	local_normal_vec.push_back(pair);
	
	i++;
      }

      // Otherwise it's 2-dim'l
      else {
	pair[0] = pow;
	if (off_diag == 0) 
	  pair[1] = -1;
	else 
	  pair[1] = -1;	

	i+=2;
      }
      
    }  // End of while loop.
       
  }  // End of p=2 case.
  
  
  // Return the local normal vector (of pairs).
  return local_normal_vec;
  
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the local conditions at p of an elementary form, described by the elementary form vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<long> LocalConditions::_Local_Conditions_for_elementary_form_at_p(const vector<long> & elementary_form_vector, const mpz_class & p) const {

  // For now, we'll use the local densities as before...
  // ---------------------------------------------------
  Matrix_mpz local_normal;
  vector<long> elementary_form_vec = elementary_form_vector;


  // Sanity Checks
  assert(p > 1);
  if (p>2)
    for(long i=0; i<elementary_form_vec.size(); i++)
      assert((elementary_form_vec[i] == 1) || (elementary_form_vec[i] == -1));
  else	     
    for(long i=0; i<elementary_form_vec.size(); i++) {
      assert(elementary_form_vec[i] >= -4);
      assert(elementary_form_vec[i] <= 14);
      assert((elementary_form_vec[i] % 4 != 0) || (elementary_form_vec[i] == -4));  
    }



  // First reorder the vector so its p-divisibility is increasing (possibly needed for the local densities routine, so we'll be safe for now).
  bool done_flag = false;
  while (done_flag == false) {
    for(long i=0; i+1 < elementary_form_vec.size(); i++)
      if ((elementary_form_vec[i] % p == 0) && (elementary_form_vec[i+1] % p != 0)) {
	long tmp = elementary_form_vec[i];
	elementary_form_vec[i] = elementary_form_vec[i+1];
	elementary_form_vec[i+1] = tmp;
	done_flag = true;
      }
  }



  // Make the (local normalized) matrix of the elementary form
  if (p>2) {
    long n = elementary_form_vec.size();
    local_normal.SetDims(n,n);
    mpz_class s = NonResidue(p);
    for(long i=0; i<n; i++)
      if (elementary_form_vec[i] >= 0) 
	local_normal(i+1,i+1) = elementary_form_vec[i];
      else 
	local_normal(i+1,i+1) = -s * elementary_form_vec[i];
  }
  else {       //  --------- This is for p=2 -------------
        
    // Compute the size of the matrix
    long n = elementary_form_vec.size();
    for(long i=0; i<elementary_form_vec.size(); i++)
      if (elementary_form_vec[i] < 0) 
	n++;
    
    // Make the new (local normalized) matrix for p=2
    local_normal.SetDims(n,n);
    long index = 1;
    for(long i=0; i<elementary_form_vec.size(); i++)

      // It's 1-dim'l
      if (elementary_form_vec[i] >= 0) {
	local_normal(index, index) = elementary_form_vec[i];
	index++;
      }
    
      // It's 2-dim'l
      else {

	// If it's xy or 2xy
	if (elementary_form_vec[i] >= -2) {
	  local_normal(index, index) = 0;
	  local_normal(index, index+1) = 1;
	  local_normal(index+1, index+1) = 0;
	}	  

	// If it's x^2 + xy + y^2 or 2(x^2 + xy + y^2)
	else {
	  local_normal(index, index) = 1;
	  local_normal(index, index+1) = 1;
	  local_normal(index+1, index+1) = 1;
	}

	// Finally, add the factor of 2 if needed
	if (elementary_form_vec[i] % 2 == 0) {
	  local_normal(index, index) *= 2;
	  local_normal(index, index+1) *= 2;
	  local_normal(index+1, index+1) *= 2;
	}
	
	// Increment the index
	index += 2;

      }
    
  }
  

  // Find the level at p
  mpz_class N = 1;
  if (p > 2) { 
    for(long i=0; i<elementary_form_vec.size(); i++)
      if ((N == 1) && (elementary_form_vec[i] % p == 0))  // Make the level p if we have any p-power (which must be p^1 since it's elementary)
	N = p;
  }
  else {   // For p = 2
    for(long i=0; i<elementary_form_vec.size(); i++) {
      if ((N < 8) && (elementary_form_vec[i] % 2 == 0) && (elementary_form_vec[i] > 0))
	N = 8;
      else if ((N < 4) && (elementary_form_vec[i] % 2 != 0) && (elementary_form_vec[i] > 0))
	N = 4;
      else if ((N < 2) && (elementary_form_vec[i] % 2 == 0) && (elementary_form_vec[i] < 0))
	N = 2;
    }
  }
  

  // Find the local conditions
  return Make_Local_Repn_Vector_At_Prime(local_normal, N, p);
  



  /*
  // Make the local conditions vector
  vector<long> local_conditions(9,-99);
  local_conditions[0] = p.get_ui();

  // Return the local conditions
  return local_conditions;
  */

}



//////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the local conditions at p of an elementary form, described by the elementary form vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////
vector<long> LocalConditions::_Local_Conditions_Q_at_p(const Matrix_mpz & QQ, const mpz_class & p) const {

  // Find the local normal form Q of QQ at p
  Matrix_mpz Q = QQ.GetLocalNormal(p);
  
  // Convert the local normal form into a vector of pairs
  vector< vector<long> > pair_vector = _Parse_Local_Normal_Form_at_p(Q,p);
  long pair_index = 1;


  // Find the overall factor of p ^ nu in Q
  unsigned long nu = max(Valuation(Q(1,1), p), Valuation(Q(1,2), p));      // Check that this works for xy at p=2

  // Declare some useful vectors
  vector<long> universal_conditions(9,0);
  universal_conditions[0] = p.get_ui();
  bool universal_flag = false;
  vector<long> local_conditions(9,0);
  local_conditions[0] = p.get_ui();
  vector<long> tmp_local_conditions(9,0);

  // Find the local conditions for Q1
  vector<long> elem_form_vec(1, pair_vector[pair_index][1]);  // Make the elementary form vector for Q1
  tmp_local_conditions = _Local_Conditions_for_elementary_form_at_p(elem_form_vec, p);
  if (tmp_local_conditions == universal_conditions)
    universal_flag == true;
  
  // Multiply them by p^(nu1)
  tmp_local_conditions = _Multiply_Local_Conditions_by_power_at_p(tmp_local_conditions, pair_vector[pair_index][0], p);
  pair_index++;
  

  // Loop if we have more elementary local forms and the (p^e)Q1 + Q2 we just added was not universal
  while ((pair_index < pair_vector.size()) && (universal_flag == false)) { 
    long nu_index = pair_vector[pair_index][0];

    // Make the elementary form vector corresponding to the (p^*)Q1 + Q2 sublattice
    elem_form_vec.clear();
    for(long i=0; i<=pair_index; i++)                 // Make the basic types
      elem_form_vec.push_back(pair_vector[i][1]);
    for(long i=0; i<pair_index; i++)                  // Multiply by p where necessary (the p-power has a different parity)
      if ((nu_index - pair_vector[i][0]) % 2 != 0)
	elem_form_vec[i] *= p.get_si();
    
    // Find the elementary conditions for (p^e)Q1 + Q2 and multiply them by p^(nu2)
    tmp_local_conditions = _Local_Conditions_for_elementary_form_at_p(elem_form_vec, p);
    if (tmp_local_conditions == universal_conditions)
      universal_flag == true;
    tmp_local_conditions = _Multiply_Local_Conditions_by_power_at_p(tmp_local_conditions, nu_index, p);
    
    // Add it to the cumulative local conditions (i.e. take the union)
    for(long i=1; i<=8; i++)
      if (local_conditions[i] < 0)
	local_conditions[i] = tmp_local_conditions[i];
      else
	local_conditions[i] = min(local_conditions[i],tmp_local_conditions[i]);

  }


  // Return the local conditions
  return local_conditions;
}




// -----------------------------------------------------------------------------------------------


///////////////////////////////////////////////////////////
// This returns the union of a vector of local conditions
////////////////////////////////////////////////////////////

LocalConditions Union(const vector<LocalConditions> & local_vec) {

  LocalConditions new_union;

  // Deal with the empty vector
  if (local_vec.size() == 0)
    return new_union;

  new_union = local_vec[0];

  // Deal with a vector of one element
  if (local_vec.size() == 1) 
    return new_union;

  // Make the union
  for(long i=1; i < local_vec.size(); i++)
    new_union.UnionWith(local_vec[i]);

  return new_union;

}



/////////////////////////////////////////////////////////////////////
// This returns the intersection of a vector of local conditions
/////////////////////////////////////////////////////////////////////

LocalConditions Intersection(const vector<LocalConditions> & local_vec) {

  LocalConditions new_intersection;

  // Deal with the empty vector
  if (local_vec.size() == 0)
    return new_intersection;

  new_intersection = local_vec[0];

  // Deal with a vector of one element
  if (local_vec.size() == 1) 
    return new_intersection;

  // Make the intersection
  for(long i=1; i < local_vec.size(); i++)
    new_intersection.IntersectWith(local_vec[i]);

  return new_intersection;

}




// ================================ Strict I/O Front-end Routines ==========================================


///////////////////////////////////////////////
// Writes the local conditions to the ostream
///////////////////////////////////////////////
void FileOutput(const LocalConditions & elt, ostream & out) {

  elt._FileOutput(out);

}


////////////////////////////////////////////////
// Reads the local conditions from the istream
////////////////////////////////////////////////
void FileInput(LocalConditions & elt, istream & in) {

  elt._FileInput(in);

}






