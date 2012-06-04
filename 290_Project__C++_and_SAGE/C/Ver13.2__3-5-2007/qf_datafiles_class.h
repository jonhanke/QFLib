

#if !defined(QF_DATAFILES_CLASS_H)
#define QF_DATAFILES_CLASS_H



class QF_Datafiles {

public:
  QF_Datafiles(const string & projectname = "Temp_QF_Project__Default");  // Constructor

  // Flag routines ... add these later if we make the flags private... =)
  

  // Directory Clearing routines
  void Clear_Eis();    // Clears out all of the Eisenstein series data files
  void Clear_Theta();    // Clears out all of the theta series data files
  void Clear_Boolean_Theta();    // Clears out all of the boolean theta series data files
  void Clear_Approx_Boolean_Theta();    // Clears out all of the approximate boolean theta series data files
  void Clear_Theta_All();  // Clears all types of theta series data files

  void Clear_Eis_Bounds();
  void Clear_Cusp_Bounds();
  void Clear_F4_Bounds();
  void Clear_Bounds_All();

  void Clear_Eligible_Primes();
  void Clear_Local_Covers();
  void Clear_Exceptions();
  void Clear_Overflows();
  void Clear_Ternary_Exceptions();

  void Clear_All();    // Clears out all of the old datafiles
  void Delete_All_Folders();


  // Public Variables:
  // -----------------

  // Directories
  string Project_Dir;   // Directory where all other files and directories are stored

  string Eis_Dir;  // Directory where the Eisenstein Series are stored
  string Theta_Dir;  // Directory where the theta functions are stored
  string Boolean_Theta_Dir; // Directory where the boolean theta functions are stored (these are not approximate!)
  string Approx_Boolean_Theta_Dir; // Directory where the *approximate* boolean theta functions are stored 

  string Eis_Lower_Bound_Dir;  // Directory where the Eisenstein coefficient lower bounds are stored
  string Cusp_Upper_Bound_Dir;  // Directory where the cuspidal coefficient upper bounds are stored
  string F4_Bound_Dir;  // Directory where the combined F4 bounds are stored

  string Eligible_Prime_Dir;  // Directory where the eligible primes are stored
  string Local_Cover_Dir;     // Directory where the (minimal) local covers are stored
  string Exceptions_Dir;     // Directory where the exceptional numbers are stored (for 4-variable forms)
  string Overflows_Dir;     // Directory where the overflowed numbers are stored (for 4-variable forms)
  string Ternary_Exceptions_Dir;     // Directory where the ternary exceptional numbers are stored


  // Flags
  bool Keep_Eligible_Primes;  // Sets whether we keep the list of eligible primes in a file
  bool Keep_Local_Covers;     // Sets whether we keep the minimal local covers in a file
  bool Keep_Exceptions;       // Sets whether we keep the exceptions in a file

  
private:

  // Variables:
  // ----------


};



// -------------------------------------------------------------------------------------


#include "qf_datafiles_class.cc"


#endif






