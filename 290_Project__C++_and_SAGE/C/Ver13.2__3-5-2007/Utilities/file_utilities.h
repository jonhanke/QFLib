

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



// ---------------------------

#if !defined(FILE_UTILITIES_H)
#define FILE_UTILITIES_H

#include "file_utilities.cc"

#endif
