

#if !defined(QF_PROJECT_CLASS_H)
#define QF_PROJECT_CLASS_H


#include "qf_datafiles_class.h"



class QF_Project : public QF_Datafiles {     // This inherits the QF_Datafiles class as well! =)

public:
  // Constructor
  static set<long> EMPTY_SET;
  QF_Project(const string & projectname = "Temp_QF_Project", const string & filename = "",
	     const string & cusp_const_dir = "", const string & cusp_const_prefix = "", const set<long> & excluded_forms = EMPTY_SET);  

  // File Operations
  void Read_Approx_Forms_and_Cusp_Bounds(const string & filename);                               
  void Read_Exact_Forms_and_Cusp_Bounds(const string & filename, const string & cusp_dir, const string & cusp_prefix);                               
  ostream & Print(ostream & out) const;   


  // Project Computations
  void FindExceptions(const string & type = "minimal", 
		      const long & first = 1, const long & last = 0);  // Finds all exceptions for the current project

  void PrintExceptions(const long & first = 1, const long & last = 0);



  // Public Variables:
  // -----------------

  // Datafiles
  QF_Datafiles Datafolders;

  // Project Lists
  set<long> Excluded_forms; 
  vector<Matrix_mpz> Form_List;
  vector<mpz_class> Level_List;  
  vector<double> Cusp_Constant_List;  // *** THIS SHOULD BE REMOVED WHEN WE HAVE MAGMA/PARI INTERFACE! =) ***
                                      //     Also the read routine will need to be changed...

  vector< vector<mpz_class> > Exception_Lists;  // This stores all of the exceptions for the project.
  set<mpz_class> Project_exception_set;         // Running list of exceptions for the project.

  
private:

  // Variables:
  // ----------


};




// Here are some extras:
// ---------------------
ostream & operator<<(ostream & out, const QF_Project & local);



// -------------------------------------------------------------------------------------


// These may or may not be necessary (silly output routines...)    <--- We should check this later...
#include "output.h"

#include "representability.h"

#include "qf_project_class.cc"


#endif






