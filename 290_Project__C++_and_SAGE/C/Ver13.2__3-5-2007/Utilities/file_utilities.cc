

/////////////////////////////
// Checks if a file exists //
/////////////////////////////

bool FileExists(char* filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool FileExists(const char* & filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool FileExists(string filename) {

  ifstream file;
  file.open(filename.c_str());

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}




//////////////////////////////////
// Checks if a directory exists //   <== NOTE: This is a copy of FileExists()! 
//////////////////////////////////

bool DirectoryExists(char* filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool DirectoryExists(const char* & filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool DirectoryExists(string filename) {

  ifstream file;
  file.open(filename.c_str());

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}




///////////////////////////////////////////////////
// Returns the absolute path to a file/directory //
///////////////////////////////////////////////////

string GetAbsolutePath(string filename) {

  // Check if the given path is already absolute
  if (filename[0] == '/') {
    //    string newfilename(filename);
    //    return newfilename;
    
    return string(filename);
  }

  // If not, then either find the path to "~/" or to "./"
  string newfilename(filename);
  string filepath;
  if ((filename[0] == '~') && (filename[1] == '/')) {
    filepath = string(getenv("HOME"));           // Find the home directory path
    newfilename.erase(0,2);                      // Forget the first two characters
  }
  else if ((filename[0] == '.') && (filename[1] == '/')) {
    filepath = string(getenv("PWD"));            // Find the present working directory path
  }
  else {
    filepath = string(getenv("PWD"));            // Find the present working directory path
  }


  // Concatenate this with the given filename
  string absolutepath;
  absolutepath = filepath + "/" + newfilename;   

  // Return the absolute path
  return absolutepath;

}


/////////////////////////////////////////////////////////////////////////////////////
// Returns the absolute path of the project directory, and creates it if necessary //
/////////////////////////////////////////////////////////////////////////////////////
string GetAbsoluteProjectDirPath() {

  // Try to get the path from the shell variable QFLIB_local_project_data_dir 
  string qflib_project_data_dir = string("");
  if(std::getenv("QFLIB_local_project_data_dir") != NULL) {
    qflib_project_data_dir = string(std::getenv("QFLIB_local_project_data_dir")) + string("/");
  }

  // Set a default path if no shell variable is found
  string default_string("/tmp/QFLIB_Project_Data/");
  if (qflib_project_data_dir.empty()) {
    cout << "Oops!  There is no QFLIB_local_project_data_dir shell variable to use -- trying the default setting" << endl;
    cout << "    Using QFLIB_local_project_data_dir = " << default_string << endl;
    qflib_project_data_dir = default_string;
  }

  // Check that this data directory exists, and create it if it doesn't exist
  if (not DirectoryExists(qflib_project_data_dir)) {
    string make_qflib_project_data_dir_string = string("mkdir -p ") + qflib_project_data_dir;
    system(make_qflib_project_data_dir_string.c_str());
  }
  
  // Return the absolute project path
  return GetAbsolutePath(qflib_project_data_dir);

}



///////////////////////////////////////////////////////////////////
// Returns the given or default path of the remote MAGMA program //
///////////////////////////////////////////////////////////////////
string GetRemoteMAGMAPath() {

  // Try to get the path from the shell variable QFLIB_remote_MAGMA_path
  string remote_magma_path = string("");
  if(std::getenv("QFLIB_remote_MAGMA_path") != NULL) {
    remote_magma_path = string(std::getenv("QFLIB_remote_MAGMA_path"));
  }

  // Set a default path if no shell variable is found
  string default_string("/usr/local/bin/magma");
  if (remote_magma_path.empty()) {
    cout << "Oops!  There is no QFLIB_remote_MAGMA_path shell variable to use -- trying the default setting" << endl;
    cout << "    Using QFLIB_remote_MAGMA_path = " << default_string << endl;
    remote_magma_path = default_string;
  }

  // Return the remote magma path
  return remote_magma_path;
  
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the given or default path of the remote temporary directory for cuspidal computations //
///////////////////////////////////////////////////////////////////////////////////////////////////
string GetRemoteTemporaryComputationsDirectory() {

  // Try to get the path from the shell variable QFLIB_remote_temp_dir
  string remote_temporary_dir_path = string("");
  if(std::getenv("QFLIB_remote_temp_dir") != NULL) {
    remote_temporary_dir_path = string(std::getenv("QFLIB_remote_temp_dir"));
  }

  // Set a default path if no shell variable is found
  string default_string("/tmp/QFLIB_Temporary_Computations/");
  if (remote_temporary_dir_path.empty()) {
    cout << "Oops!  There is no QFLIB_remote_temp_dir shell variable to use -- trying the default setting" << endl;
    cout << "    Using QFLIB_remote_temp_dir = " << default_string << endl;
    remote_temporary_dir_path = default_string;
  }
  
  // Return the remote temporary computations directory path
  return remote_temporary_dir_path;

}



////////////////////////////////////////////////////////////////////////////////////////
// Returns the user@machine string used for remote computations of cuspidal constants //
////////////////////////////////////////////////////////////////////////////////////////
string GetRemoteUserAtMachine() {

  // Try to get the path from the shell variable QFLIB_remote_user_at_machine 
  string remote_user_at_machine = string("");
  if(std::getenv("QFLIB_remote_user_at_machine") != NULL) {
    remote_user_at_machine = string(std::getenv("QFLIB_remote_user_at_machine"));
  }

  // Fail if no shell variable is found
  if (remote_user_at_machine.empty()) {
    cout << "Oops!  There is no QFLIB_remote_user_at_machine shell variable to use -- set this!" << endl;
    assert(0==1);
  }

  // Return the remote user@machine string
  return remote_user_at_machine;
  
}





/////////////////////////////////////////////////////////////////////////
// Returns the absolute path to the directory holding the primes files //
/////////////////////////////////////////////////////////////////////////
string GetPrimesDirectory() {

  // Try to get the path from the shell variable QFLIB_Project_Data_Dir 
  string primes_dir_path = string("");
  if(std::getenv("QFLIB_local_primes_dir") != NULL) {
    primes_dir_path = string(std::getenv("QFLIB_local_primes_dir"));
  }

  cout << "Inside GetPrimesDirectory() -- step 1" << endl;  

  // Set a default path if no shell variable is found
  string default_string("../../Primes/");
  cout << "Inside GetPrimesDirectory() -- step 1a" << endl;  
  if (primes_dir_path.empty()) {
    cout << "Inside GetPrimesDirectory() -- step 1b" << endl;  
    cout << "Oops!  There is no QFLIB_local_primes_dir shell variable to use -- trying the default setting" << endl;
    cout << "    Using QFLIB_local_primes_dir = " << default_string << endl;
    primes_dir_path = default_string;
  }

  cout << "Inside GetPrimesDirectory() -- step 2" << endl;  
  
  // Return the remote temporary computations directory path
  return GetAbsolutePath(primes_dir_path);

}

