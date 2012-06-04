


// Constructor -- requires a bounded number of eligible prime factors
LocalCheckingStats::LocalCheckingStats(const long & num_of_prime_factors){
  
  MAX_STATS_DEPTH = 100;       // This keeps stats for at most 100 depths (from 0 to 99).

  long buffered_num_of_primes = num_of_prime_factors + 1;  // This lets the indexing agree with the # of prime factors

  // Sanity Check
  assert(num_of_prime_factors > 0);
  assert(num_of_prime_factors < 30);
  
  // Clear the existing vectors
  total_eligible.clear();
  total_represented.clear();
  total_missed.clear();
  representation_depths.clear();  

  // Resize and initialize them to zero
  total_eligible.resize(buffered_num_of_primes, 0);
  total_represented.resize(buffered_num_of_primes, 0);
  total_missed.resize(buffered_num_of_primes, 0);
  vector<mpz_class> empty_vec(MAX_STATS_DEPTH + 1, 0);     
  representation_depths.resize(buffered_num_of_primes, empty_vec);   // This is initialized to use the default depth 

}



// Empty Destructor 
LocalCheckingStats::~LocalCheckingStats(){
  // Do nothing...
}



// ----------------------------------------------------------------------------------------------

/*
// Prints an mpz_class as a right formated object, in a certain amount of space! 
ostream & Rprint LocalCheckingStats::Print(ostream & out, mpz_class m) const { 
}
*/


////////////////////////////////
// Prints the local conditions
///////////////////////////////
ostream & LocalCheckingStats::Print(ostream & out) const {    

  // Show total missed/represented statistics
  cout << "   # of primes            # eligible            # represented        # missed       % missed " << endl;
  cout << "   -----------            ----------            -------------        --------       -------- " << endl;
  for (long num = max(0,total_eligible.size() - 1); num >= 1; num--) {
    cout << "       " << setw(2) << num << "            " 
	 << setw(12) << total_eligible[num] << "            " 
	 << setw(12) << total_represented[num] << "           " 
	 << setw(6) << total_missed[num] << "            "; 
    if (total_eligible[num] > 0) {
      cout.setf(ios_base::fixed, ios_base::floatfield);
      cout.precision(1);
      cout << (100.0 * total_missed[num].get_d() / total_eligible[num].get_d()) << "%   "; 
      cout.setf(ios_base::fmtflags(0), ios_base::floatfield);
      cout.precision(6);
    }
    else 
      cout << "XXX  ";
    cout << endl;
  }
  cout << endl << endl;    

  
  // Compute the maximum depths for each number of prime factors
  vector<long> max_depth_vec(total_eligible.size());
  for (long num = 1; num < total_eligible.size(); num++) 
    for(long depth_index = representation_depths[num].size() - 1; depth_index >= 0; depth_index--) 
      if ((max_depth_vec[num] == 0) && (representation_depths[num][depth_index] > 0))
	max_depth_vec[num] = depth_index;
  

  // Precompute the size of the largest number appearing for each depth
  vector<long> depth_count__string_size(MAX_STATS_DEPTH, 0);  // Be sure to initialize these values to zero! =)
  for(long i=0; i < MAX_STATS_DEPTH; i++)
    for(long j=0; j < total_eligible.size(); j++) {
      
      // Find the number of digits in a number
      mpz_class num(representation_depths[j][i]);
      long num_of_digits=0;
      mpz_class tmp_num(1);
      while (tmp_num <= num) {
	tmp_num = 10 * tmp_num;
	num_of_digits++;
      }

      // Take the maximum number of digits for each depth
      depth_count__string_size[i] = max(num_of_digits, depth_count__string_size[i]);
    }  
  


  // Show depth % distibution statistics
  cout << " Depth distribution percentages (starting at depth 0): " << endl;
  cout << " ----------------------------------------------------- " << endl;
  for (long num = 1; num < total_eligible.size(); num++) {
    long depth = 0;
    cout << " # of primes = " << setw(2) << num << ":   ";
    while (depth <= max_depth_vec[num]) {   
      if (total_represented[num] > 0) {
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout.precision(1);
	cout << setw(5) << (100.0 * representation_depths[num][depth].get_d() / total_represented[num].get_d()) << "%   ";
	cout.setf(ios_base::fmtflags(0), ios_base::floatfield);
	cout.precision(6);
      }
      else 
	cout << "XXX  ";
      depth++;
    }
    cout << endl;
  }
  cout << endl << endl;  


  // Show depth count statistics
  cout << " Depth distributions counts (starting at depth 0): " << endl;
  cout << " ------------------------------------------------- " << endl;
  for (long num = 1; num < total_eligible.size(); num++) {
    long depth = 0;
    cout << " # of primes = " << setw(2) << num << ":   ";
    while (depth <= max_depth_vec[num]) {   
      cout << setw(depth_count__string_size[depth] + 2) << representation_depths[num][depth] << "   "; 
      depth++;
    }
    cout << endl;
  }
  cout << endl << endl;  


  // Show the maximum depths necessary for each number of prime factors
  cout << " Maximum necessary depths (=< 100) for represented eligible numbers with a given number of prime factors: " << endl;
  cout << " -------------------------------------------------------------------------------------------------------- " << endl;
  for (long num = 1; num < total_eligible.size(); num++) 
        cout << " # of primes = " << setw(2) << num << ":     "
	 << " maximum depth needed = " << setw(3) << max_depth_vec[num] << endl; 

  
  cout << endl << endl;  

  /*
  // Show the maximum depths necessary for each number of prime factors
  cout << " Maximum necessary depths (=< 100) for missed eligible numbers with a given number of prime factors: " << endl;
  cout << " --------------------------------------------------------------------------------------------------- " << endl;
  for (long num = 0; num < total_eligible.size(); num++) {

    long max_depth = 0;
    for(long depth_index = representation_depths[num].size(); depth_index>0; depth_index--) {

      // Find the maximum depth
      if ((max_depth == 0) && (representation_depths[num][depth_index-1] > 0))
	max_depth = representation_depths[num][depth_index-1];

      // Print them
      cout << " # of primes = " << setw(2) << num << ":     "
	   << " maximum depth needed = " << setw(3) << max_depth << endl; 

    }

  }
  cout << endl << endl;  
  */


  return out;
}


////////////////////////////////////////////////////
// This is the overloaded (global) version for <<
////////////////////////////////////////////////////
ostream & operator<<(ostream & out, const LocalCheckingStats & local) {    
  return local.Print(out);
}








///////////////////////////////////////////////
// Writes the local conditions to the ostream
///////////////////////////////////////////////
void LocalCheckingStats::_FileOutput(ostream & out) const {

  /*

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

  */

}



////////////////////////////////////////////////
// Reads the local conditions from the istream
////////////////////////////////////////////////
void LocalCheckingStats::_FileInput(istream & in) {

  /*

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

  */

}


// ================================ Strict I/O Front-end Routines ==========================================


///////////////////////////////////////////////
// Writes the local conditions to the ostream
///////////////////////////////////////////////
void FileOutput(const LocalCheckingStats & elt, ostream & out) {

  elt._FileOutput(out);

}


////////////////////////////////////////////////
// Reads the local conditions from the istream
////////////////////////////////////////////////
void FileInput(LocalCheckingStats & elt, istream & in) {

  elt._FileInput(in);

}






