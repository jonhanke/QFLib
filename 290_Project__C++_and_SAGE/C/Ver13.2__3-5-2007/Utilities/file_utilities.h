

// Check if a file exists.
bool FileExists(char* filename);
bool FileExists(const char* & filename);
bool FileExists(string filename);
//bool FileExists(const string & filename);

// Check if a directory exists. 
bool DirectoryExists(char* filename);
bool DirectoryExists(const char* & filename);
bool DirectoryExists(string filename);
//bool DirectoryExists(const string & filename);

// Get the absolute path to a file/directory
string GetAbsolutePath(string filename);

// Get the absolute path to the QFLIB_Project_Data directory
string GetAbsoluteProjectDirPath();

// Get the path of the remote MAGMA executable
string GetRemoteMAGMAPath();

// Get the path of the remote temporary directory
string GetRemoteTemporaryComputationsDirectory();

// Get the user@machine for remote cuspidal computations
string GetRemoteUserAtMachine();

// Get the absolute path to the local Primes directory
string GetPrimesDirectory();


// ---------------------------

#if !defined(FILE_UTILITIES_H)
#define FILE_UTILITIES_H

#include "file_utilities.cc"

#endif
