
// Constructor to make a project directory system (if it doesn't exist) 
QF_Datafiles::QF_Datafiles(const string & projectname) {

  extern const char ABSOLUTE_PROJECT_PATH[25];

  // Make the project directory name
  char absolute_project_path[300];
  sprintf(absolute_project_path, "%s%s/", GetAbsolutePath(ABSOLUTE_PROJECT_PATH).c_str(), projectname.c_str());     // This also appends the trailing '/'. =)


  // Check if "projectname" exists in "~/QF_Project_Data/", else create a project
  if (DirectoryExists(absolute_project_path) == false) {
    char command_line[400];
    
    // Make a new project directory
    sprintf(command_line, "mkdir %s", absolute_project_path);
    system(command_line);

    // Make a new Eis_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Eis_Dir/");
    system(command_line);

    // Make a new Theta_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Theta_Dir/");
    system(command_line);

    // Make a new Boolean_Theta_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Boolean_Theta_Dir/");
    system(command_line);

    // Make a new Approx_Boolean_Theta_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Approx_Boolean_Theta_Dir/");
    system(command_line);

    // Make a new Eis_Lower_Bound_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Eis_Lower_Bound_Dir/");
    system(command_line);

    // Make a new Cusp_Upper_Bound_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Cusp_Upper_Bound_Dir/");
    system(command_line);

    // Make a new F4_Bound_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "F4_Bound_Dir/");
    system(command_line);

    // Make a new Eligible_Primes_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Eligible_Primes_Dir/");
    system(command_line);

    // Make a new Local_Cover_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Local_Cover_Dir/");
    system(command_line);

    // Make a new Exceptions_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Exceptions_Dir/");
    system(command_line);

    // Make a new Overflows_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Overflows_Dir/");
    system(command_line);

    // Make a new Ternary_Exceptions_Dir
    sprintf(command_line, "mkdir %s%s", absolute_project_path, "Ternary_Exceptions_Dir/");
    system(command_line);

  }


  // Set the absolute directory names
  Project_Dir = absolute_project_path;

  Eis_Dir = Project_Dir + "Eis_Dir/";
  Theta_Dir = Project_Dir + "Theta_Dir/";
  Boolean_Theta_Dir = Project_Dir + "Boolean_Theta_Dir/"; 
  Approx_Boolean_Theta_Dir = Project_Dir + "Approx_Boolean_Theta_Dir/";

  Eis_Lower_Bound_Dir = Project_Dir + "Eis_Lower_Bound_Dir/";
  Cusp_Upper_Bound_Dir = Project_Dir + "Cusp_Upper_Bound_Dir/";
  F4_Bound_Dir = Project_Dir + "F4_Bound_Dir/";

  Eligible_Prime_Dir = Project_Dir + "Eligible_Prime_Dir/";
  Local_Cover_Dir = Project_Dir + "Local_Cover_Dir/";
  Exceptions_Dir = Project_Dir + "Exceptions_Dir/";
  Overflows_Dir = Project_Dir + "Overflows_Dir/";
  Ternary_Exceptions_Dir = Project_Dir + "Ternary_Exceptions_Dir/";


  // Set flag to their default values
  Keep_Eligible_Primes = false;  // Don't keep the list of eligible primes in a file
  Keep_Local_Covers = true;      // Keep the minimal local covers in a file
  Keep_Exceptions = true;        // Keep the exceptions in a file

}



// ==================================== Directory Clearing Routines =================================================


void QF_Datafiles::Clear_Eis() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Eis_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Theta() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Theta_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Boolean_Theta() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Boolean_Theta_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Approx_Boolean_Theta() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Approx_Boolean_Theta_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Theta_All() {
  Clear_Theta();
  Clear_Boolean_Theta();
  Clear_Approx_Boolean_Theta();
}


// -----------------------------------------


void QF_Datafiles::Clear_Eis_Bounds() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Eis_Lower_Bound_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Cusp_Bounds() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Cusp_Upper_Bound_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_F4_Bounds() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", F4_Bound_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Bounds_All() {
  Clear_Eis_Bounds();
  Clear_Cusp_Bounds();
  Clear_F4_Bounds();
}


// -----------------------------------------


void QF_Datafiles::Clear_Eligible_Primes() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Eligible_Prime_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Local_Covers() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Local_Cover_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Exceptions() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Exceptions_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Overflows() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Overflows_Dir.c_str());
  system(command_line);
}


void QF_Datafiles::Clear_Ternary_Exceptions() {
  char command_line[400];
  sprintf(command_line, "rm -f  %s*", Ternary_Exceptions_Dir.c_str());
  system(command_line);
}


// -----------------------------------------


void QF_Datafiles::Clear_All() {
  Clear_Theta_All();
  Clear_Bounds_All();

  Clear_Eligible_Primes();
  Clear_Local_Covers();
  Clear_Exceptions();
  Clear_Overflows();
  Clear_Ternary_Exceptions();
}


void QF_Datafiles::Delete_All_Folders() {
  char command_line[400];
  sprintf(command_line, "rm -rf  %s", Project_Dir.c_str());
  system(command_line);
}


// ====================================================================================


