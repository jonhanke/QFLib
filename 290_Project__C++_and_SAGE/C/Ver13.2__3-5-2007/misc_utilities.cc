
#include <gmp.h>
#include <gmpxx.h>

#include "LCM.h"


/*
//  <<====  THIS IS ALREADY DEFINED IN THE EARLIER LOCAL DENSITY PROJECT...

// LCM of two mpz_classes
mpz_class LCM(mpz_class a, mpz_class b) {

  mpz_t l;
  mpz_init(l);
  mpz_lcm(l, a.get_mpz_t(), b.get_mpz_t());

  // Move the LCM
  mpz_class lcm;
  lcm = mpz_class(l);

  // Clear the mpz_t variable
  mpz_clear(l);

  return lcm;
}
*/


// -------------------------------------------------------


// Print the current date and time 
void PrintTime() {
  time_t rawtime;
  struct tm * timeinfo;
  
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  cout << " " << asctime(timeinfo) << endl;
  //  cout << "Current date and time are: " << asctime(timeinfo) << endl;
}


// Return the current date/time
string Time() {
  time_t rawtime;
  struct tm * timeinfo;
  
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  return asctime(timeinfo);
}




/*
// Define the << operator for the vector<T> type
template <class T>
ostream & operator<<(ostream & out, const vector<T> & v) {    
  out << "[ ";
  for(long i=0; i < v.size(); i++)
    out << v[i] << " ";
  out << " ]";
  return out;
}
*/



// --------------------------------------------------------------




// Finds the level of a ternary form Q
long QF_Ternary_Level(long Q[6]) {

  // Make a 3 x 3 matrix
  mpq_class M[3][3];
  M[0][0] = mpq_class(Q[0]);
  M[0][1] = mpq_class(Q[1]) / 2;
  M[1][0] = mpq_class(Q[1]) / 2;
  M[0][2] = mpq_class(Q[2]) / 2;
  M[2][0] = mpq_class(Q[2]) / 2;
  M[1][1] = mpq_class(Q[3]);
  M[1][2] = mpq_class(Q[4]) / 2;
  M[2][1] = mpq_class(Q[4]) / 2;
  M[2][2] = mpq_class(Q[5]);

  // Double the matrix
  for (long i=0; i<3; i++)
    for (long j=0; j<3; j++)
      M[i][j] = 2 * M[i][j];

  // Find the determinant of the matrix
  mpq_class M_det;
  M_det = M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2])
    - M[0][1] * (M[1][0] * M[2][2] - M[2][0] * M[1][2]) 
    + M[0][2] * (M[1][0] * M[2][1] - M[2][0] * M[1][1]); 

  // Find the adjoint of the matrix
  mpq_class Adj[3][3];
  Adj[0][0] = (M[1][1] * M[2][2] - M[2][1] * M[1][2]);
  Adj[0][1] = -(M[0][1] * M[2][2] - M[2][1] * M[0][2]);
  Adj[1][0] = -(M[1][0] * M[2][2] - M[1][2] * M[2][0]);
  Adj[0][2] = (M[0][1] * M[1][2] - M[1][1] * M[0][2]);
  Adj[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  Adj[1][1] = (M[0][0] * M[2][2] - M[2][0] * M[0][2]);
  Adj[1][2] = -(M[0][0] * M[1][2] - M[1][0] * M[0][2]);
  Adj[2][1] = -(M[0][0] * M[2][1] - M[0][1] * M[2][0]);
  Adj[2][2] = (M[0][0] * M[1][1] - M[1][0] * M[0][1]);
  
  // Find the inverse of the matrix
  mpq_class Inv[3][3];
  for (long i=0; i<3; i++)
    for (long j=0; j<3; j++) {
      Inv[i][j] = Adj[i][j] / M_det;
      Inv[i][j].canonicalize();
    }
  
  // Find the LCM of the denominators of the inverse
  mpz_class Max_denom = 1;
  mpz_class temp_denom;
  for (long i=0; i<3; i++)
    for (long j=0; j<3; j++) {
      temp_denom = abs(Inv[i][j].get_den());
      Max_denom = LCM(Max_denom, temp_denom);
    }

  Max_denom = 2 * Max_denom;  // Throw in an extra 2 to get even diagonals...

  /*
  cout << "Max_denom = " << Max_denom << endl;  
  */

  // Make the global level variable
  mpz_class lvl;
  
  // Run through all divisors of Max_denom to determine the level
  for (mpz_class m=1; m<=Max_denom; m++) {
    bool is_valid_level = true;
    
    // Check if m | Max_denom
    if (Max_denom % m == 0) {

      /*
      cout << m << " is a divisor of " << Max_denom << endl;
      */      

      // Check all entries are integral
      for (long i=0; i<3; i++)
	for (long j=0; j<3; j++) {
	  mpq_class tmp;
	  tmp = mpq_class(m) * Inv[i][j];
	  tmp.canonicalize();
	  //	  cout << "tmp = " << tmp << endl;
 	  if (abs(tmp.get_den()) != 1)
	    is_valid_level = false;
	}
  
      //      cout << " integral? " << is_valid_level << endl;
    
      // Check all diagonal entries are 2 * integral
      for (long i=0; i<3; i++) {
	mpq_class tmp; 
	tmp = mpq_class(m) * Inv[i][i] / 2;
	tmp.canonicalize();
	//	cout << "diag tmp = " << tmp << endl;
	if (abs(tmp.get_den()) != 1)
	  is_valid_level = false;
      }

      //      cout << " even diagonal? " << is_valid_level << endl;
      
      // Check if the new level is valid
      if (is_valid_level == true)
	lvl = m;

    }
   
  }

  return lvl.get_ui();
  
}



/////////////////////////////////////////////////////////////
// Computes the valuation of an integer m mod some prime p //
/////////////////////////////////////////////////////////////
// Note: Assumes the given p is a prime.
 
unsigned long Valuation(const mpz_class & m, const long & p)
{
  
  // Sanity Checks
  assert( m != 0 );
  assert( p > 0 );       // Not strictly necessary, but this should be true...

  
  // Divide by p until we're done!
  unsigned long val = 0;
  mpz_class m1;
  m1 = m;
  while(m1 % p == 0) {
    val++;
    m1 = m1 / p;
  }
 
  return val;
}


// Takes a power of a long
long LongPow(long a, unsigned long b) {

  long c = 1;
  for (unsigned long i=1; i<=b; i++)
    c = c * a;

  return c;
}



// Returns a vector of prime divisors of m, in increasing order            // THERE'S AN OVERLOAD PROBLEM WITH THE GMP_CLASS_EXTRAS VERSION!!!
vector<long> PrimeDivisors1(mpz_class m) {

  // Sanity Check 
  m = abs(m);
  assert(m != 0);


  extern vector<long> Big_Prime_List;    // This is the global prime list.

  vector<long> prime_divisor_list;       // Place to put the prime divisors

  long ind=0;


  


  while((ind < Big_Prime_List.size()) && (m > 1)) {
    
    long p = Big_Prime_List[ind];

    // Check if p|m,
    if (m % p == 0) {

      // Add it to the list
      prime_divisor_list.push_back(p);

      // Remove as many powers of p as possible
      while (m % p == 0)
	m = m/p;
    }

    ind++;

  }


  // Return if we're done
  if (m==1)
    return prime_divisor_list;


  // Otherwise, there's an overflow
  else {
    cout << " PrimeDivisors ERROR: Prime list overflow for m = " << m << endl;
    assert(0==1);
  }    


}







// Make all subsets of a given set S, containing exactly num elements.
template <class T>
vector< set<T> > MakeSubsets(const set<T> & S, const long & num) {    

  // Sanity Checks
  assert(num >= 0);
  assert(num <= S.size());


  // Make the subset vector  
  vector< set<T> > subset_vec;

  
  // Deal with the empty set (in 2 cases)
  if ((S.empty() == true) || (num == 0)) {
    return subset_vec;
  }


  // Copy the set into a vector for easy indexing
  vector<T> V;
  for(typename set<T>::iterator i=S.begin(); i != S.end(); i++)
    V.push_back(*i);


  // Initialize the counting vector
  vector<long> index_vec(num,0);  
  for(long i=0; i<num; i++)
    index_vec[i] = i;
  
  long index = num-1;       // Start incrementing at the last place


  // Loop until we have all subsets
  while (index >= 0) {

    // DIAGNOSTIC
    cout << " V = " << V << endl;
    cout << " indexing_vec = " << index_vec << endl;


    // Make the associated subset
    set<T> new_subset;
    for(long j=0; j<num; j++)
      new_subset.insert(V[index_vec[j]]);

    // Add it to the list
    subset_vec.push_back(new_subset);

    // Increment
    bool inc_flag = false;      // Did we do an increment?
    while ((inc_flag == false) && (index >= 0))
      if (index_vec[index] < ((long) V.size() - num + index)) {   // This is always >= 0, and is the maximum allowed value in index_vec[index].
	index_vec[index]++;                    // Increment this index,
	for(long k=index+1; k<num; k++)        // and make all larger ones as small as possible!
	  index_vec[k] = index_vec[index] + (k - index);
	inc_flag = true;
      }
      else
	index--;

  }


  return subset_vec;
} 
