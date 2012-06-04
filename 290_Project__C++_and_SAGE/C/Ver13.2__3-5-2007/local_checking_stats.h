

#if !defined(LOCAL_CHECKING_STATS_H)
#define LOCAL_CHECKING_STATS_H


#include "output.h"  // Needed for output and strict I/O routines...


// This class tells us the checking statistics for a given local cover
class LocalCheckingStats {
  
 public:
  LocalCheckingStats(const long & num_of_prime_factors);  // Constructor
  ~LocalCheckingStats();  // Destructor

  long MAX_STATS_DEPTH;   // This is the maximum depth recorded in this class!
  
  
  ostream & Print(ostream & out) const;   // Prints the local cover checking statistics
  
  
  
  // Variables giving the stats for number checking, indexed by the number of prime factors  
  vector<mpz_class> total_eligible;
  vector<mpz_class> total_represented;
  vector<mpz_class> total_missed;
  
  // Vector giving the # of eligible numbers verified at a given depth (for each prime factor)
  vector< vector<mpz_class> > representation_depths;  
  
  

  
  // Rigid I/O routines for making persistent datafiles (called by front-end functions for consistency):
  // ---------------------------------------------------------------------------------------------------
  void _FileOutput(ostream & out) const;
  void _FileInput(istream & in);
  
};
 

ostream & operator<<(ostream & out, const LocalConditions & local);



// ------------------------------------
  
#include "local_checking_stats.cc"
  

#endif






