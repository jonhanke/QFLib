

#if !defined(DISTRIBUTED_CLASS_H)
#define DISTRIBUTED_CLASS_H


#include <vector>
#include <string>
#include <iostream>

using namespace std;



// template <class T>
class distributed_info {

public:
  distributed_info(const vector<string> & host_list);  // Constructor (given an array of hosts)
  ~distributed_info();  // Destructor
  
  string* hosts;
  long num_of_hosts;
  long* current_host_slice;     // Stores the slice currently being worked on for each host.
  char** host_assignment_time;   // Stores the time each host was assigned its slice.

  //  T current_host_info[];   // In our case, the job info will just be a long.
                               // This should come equipped with an operator<< method.



  /*
    // To Do Later:
    // -----------
    add_hosts(string hosts[]);
    remove_hosts(string hosts[]);
    deactivate_hosts(string hosts[]);
  */


private:

};



// ----------------------------------------


#include "distributed_class.cc"

#endif
