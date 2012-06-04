


// Constructor -- initializes hosts with host_list
distributed_info::distributed_info(const vector<string> & host_list) {
  
  num_of_hosts = host_list.size();
  hosts = new string[num_of_hosts];
  for(long i=0; i < num_of_hosts; i++)
    hosts[i] = host_list[i];
  current_host_slice = new long[num_of_hosts];
  host_assignment_time = new char*[num_of_hosts];

  // Write the hosts in order
  cout << "host_list.size() = "<< host_list.size() << endl;
  for(long i=0; i < num_of_hosts; i++)
    cout << "Host #" << i << " is " << hosts[i] << endl;

}


// Destructor 
distributed_info::~distributed_info() {

  delete[] current_host_slice;

}

