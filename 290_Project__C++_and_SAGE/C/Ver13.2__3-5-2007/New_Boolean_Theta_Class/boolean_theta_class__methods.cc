#include <algorithm>



// Constructor -- initializes the binary_theta with a form and a given precision
boolean_theta::boolean_theta()
{

  // Set the Precision
  _theta_precision = 0;

  // Initialize the theta function, and set the correct directories
  _Initialize();

  // Initialize the quadratic form
  Matrix_mpz Q(1,1);
  QQ = Q;

}



// Constructor -- initializes the binary_theta with a form and a given precision
boolean_theta::boolean_theta(mpz_class new_precision)
{

  // SANITY CHECKS
  assert(_theta.max_size() >= 1000000000);        // This should be a little smaller than the maximum size for the vector!
  assert(new_precision <= mpz_class("32000000000"));  // Check that we haven't exceeded  32 * 10^9  =  3.2 * 10^10


  // Set the Precision
  _theta_precision = new_precision;

  // Initialize the theta function, and set the correct directories
  _Initialize();

  // Initialize the quadratic form
  Matrix_mpz Q(1,1);
  QQ = Q;

}


// Constructor -- initializes the binary_theta with a form and a given precision
boolean_theta::boolean_theta(const Matrix_mpz & Q, mpz_class new_precision)
{

  // SANITY CHECKS
  assert(_theta.max_size() >= 1000000000);        // This should be a little smaller than the maximum size for the vector!
  assert(new_precision <= mpz_class("32000000000"));  // Check that we haven't exceeded  32 * 10^9  =  3.2 * 10^10


  // Set the Precision
  _theta_precision = new_precision;

  // Initialize the theta function, and set the correct directories
  _Initialize();

  // Initialize the quadratic form
  QQ = Q;

}


////////////////////////////////////////////////////////////////////////////////
// Initialize the theta function by allocaing space, and setting the variables
////////////////////////////////////////////////////////////////////////////////
void boolean_theta::_Initialize() {

  // Only call from constructor
  // Create and zero out the _theta vector
  _theta.erase(_theta.begin(), _theta.end());
  _theta.resize(_length(), 0);

  // Initialize the directories
  THETA_DIR = (char*)"/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/Theta_Data/";
  DIST_DIR = (char*)"/home/postdoc/jonhanke/290_Project/Check_constants/Coefficient_Checking/Eligible_Check/C/TEMP_Distributed_Dir/";  
  EXEC_DIR = (char*)"/home/postdoc/jonhanke/";

}




// TODO: Make a constructor to initialize this from a Magma file!



// Tells the highest meaningful theta coefficient
inline mpz_class boolean_theta::precision() const {
    return _theta_precision;
}


// Returns the value of the theta function at i (if it's in range)
inline bool boolean_theta::get_value(mpz_class i) const {

  // Sanity Check
  assert(i >= 0);         // Note: This is needed to ensure that the quotient and remainder are positive as well!

  // Sanity Check -- Truncation of the index to unsigned long
  if (mpz_class(i/32).fits_uint_p() == false) {
    cout << " OVERFLOW ERROR in boolean_theta::get_value(): The index " << mpz_class(i/32) << " is too large to fit in an unsigned long!" << endl;
    assert(mpz_class(i/32).fits_uint_p() == true);
  }

  // Sanity Check -- The index exceeds the precision
  if (i > (*this).precision()) {
    cout << "Error in boolean_theta::get_value(): " << endl;
    cout << "  The desired entry  i = " << i << " is greater than the precision " << precision() << endl;
    assert(i > (*this).precision());
  }		       


  // Return the boolean value
  unsigned long index = mpz_class(i/32).get_ui();
  return ((_theta[index] >> mpz_class(i % 32).get_ui()) % 2);  


  // ------------------------------------------------------------------------------------------------------------------------------
  //    mpz_fdiv_q_2exp(mpz_class(32).get_mpz_t(), i.get_mpz_t(), index);       // Compute the index, which is i / 32 (rounded down)  
  //    return ((_theta[i >> 5] >> (i % 32)) % 2);  
  // Should change this to use log_2(sizeof(long)) + 3 instead of 5 (since there are 8 bits in a byte)

}


// Sets the value of the theta function at i (if it's in range)
inline void boolean_theta::set_value(mpz_class i) {
  
  // Sanity Check
  assert(i >= 0);         // Note: This is needed to ensure that the quotient and remainder are positive as well!

  // Sanity Check -- Truncation of the index to unsigned long
  if (mpz_class(i/32).fits_uint_p() == false) {
    cout << " OVERFLOW ERROR in boolean_theta::set_value(): The index " << mpz_class(i/32) << " is too large to fit in an unsigned long!" << endl;
    assert(mpz_class(i/32).fits_uint_p() == true);
  }

  // Sanity Check -- The index exceeds the precision
  if (i > (*this).precision()) {
    cout << "Error in boolean_theta::set_value(): " << endl;
    cout << "  The desired entry  i = " << i << " is greater than the precision " << precision() << endl;
    assert(i > (*this).precision());
  }		       


  /*
  // DIAGNOSTIC
  cout << " Passed in the value " << i << endl;
  cout << " value / 32 = " << (i/32) << endl;
  cout << " Trying to set some bit in the index " << index << endl;
  cout << endl;
  */
  

  // Set the boolean value
  unsigned long index = mpz_class(i/32).get_ui();
  _theta[index] = _theta[index] | (1 << mpz_class(i % 32).get_ui());  // Set the appropriate bit 


  // -------------------------------------------------------------------------------------------------------------------------------
  //    mpz_fdiv_q_2exp(mpz_class(32).get_mpz_t(), i.get_mpz_t(), index);       // Compute the index, which is i / 32 (rounded down)
  //      _theta[i >> 5] = _theta[i >> 5] | (1 << (i % 32));  // Set the appropriate bit 
  // Should change this to use log_2(sizeof(long)) + 3 instead of 5 (since there are 8 bits in a byte)
  
}



// Clears the value of the theta function at i (if it's in range)
inline void boolean_theta::clear_value(mpz_class i) {
  
  // Sanity Check
  assert(i >= 0);         // Note: This is needed to ensure that the quotient and remainder are positive as well!

  // Sanity Check -- Truncation of the index to unsigned long
  if (mpz_class(i/32).fits_uint_p() == false) {
    cout << " OVERFLOW ERROR in boolean_theta::clear_value(): The index " << mpz_class(i/32) << " is too large to fit in an unsigned long!" << endl;
    assert(mpz_class(i/32).fits_uint_p() == true);
  }

  // Sanity Check -- The index exceeds the precision
  if (i > (*this).precision()) {
    cout << "Error in boolean_theta::clear_value(): " << endl;
    cout << "  The desired entry  i = " << i << " is greater than the precision " << precision() << endl;
    assert(i > (*this).precision());
  }		       


  // Clear the boolean value
  unsigned long index = mpz_class(i/32).get_ui();
  _theta[index] = _theta[index] & ~(1 << mpz_class(i % 32).get_ui());   // Clear the appropriate bit 

    
  //    mpz_fdiv_q_2exp(mpz_class(32).get_mpz_t(), i.get_mpz_t(), index);       // Compute the index, which is i / 32 (rounded down)
  //      _theta[i >> 5] = _theta[i >> 5] & ~(1 << (i % 32));  // Clear the appropriate bit 
  // Should change this to use log_2(sizeof(long)) + 3 instead of 5 (since there are 8 bits in a byte)
  
}


// Returns the quadratic form
Matrix_mpz boolean_theta::get_qf(){
  return QQ;
}


// Tells the number of bytes used in _theta
inline unsigned long boolean_theta::_length() const{
  return ( mpz_class((*this).precision() / 32).get_ui()) + 1;  
  // Should change this to use log_2(sizeof(long)) + 3 instead of 5 (since there are 8 bits in a byte)
}




/////////////////////////////////////////////////////////////////
// Reads the boolean ternary theta function from a binary file //
/////////////////////////////////////////////////////////////////

// Assumes a quadratic form has been set.

bool boolean_theta::read(const string & filename) {
  //ReadTernaryThetaBinary(unsigned long _theta[], long form[6], unsigned long precision) {

  // *** WARNING:  We assume that theta_bin has size precision/32 
  // ***           (which could also be coded as precision >> 5).

  // Notation:
  // ---------
  // precision of the theta function ==> ... + O(x^(precision + 1))

  // ---------------------------------------------------------------------


  // Make the new filename: 
  // ----------------------
  string newfilename;
  if (filename == "") 
    // No filename ==> Use "THETA_DIR/boolean_theta__'QQ'__'precision'.bindata"
    newfilename = string(THETA_DIR) + "boolean_theta__" + QQ.QF_String() + "__" + MakeString(precision()) + ".bindata";
  else
    newfilename = GetAbsolutePath(filename);  

  // Sanity Check that the file exists
  if ( FileExists(newfilename) == false )
    return false;

  /*
  cout << " boolean_ternary_theta.read(): Reading with the filename " << endl
       << newfilename << endl;
  */


  // -------------------------------------------------------------------

  
  // Open the file for reading and check it opened correctly
  ifstream filein;
  filein.open(newfilename.c_str(), ios::binary);


  // Flag to reflect the success of the file read
  bool file_open = true;                                 // WARNING: We should really make this flag reflect the success of the read operation too...
  if (! filein.is_open()) { 
    cout << "Error opening input file " << newfilename; 
    file_open = false; 
  }


  // Read the array
  //_theta.erase(_theta.begin(), _theta.end());
  unsigned long *buffer = new unsigned long[_length()];
  filein.read( (char*) buffer, sizeof(unsigned long) * _length() );
  _theta.assign(buffer, buffer + _length());
  delete [] buffer;
  buffer = 0;

  // Close the file
  filein.close();


  return file_open;
}






////////////////////////////////////////////////////////
// Writes the boolean theta function as a binary file //
////////////////////////////////////////////////////////

void boolean_theta::write(const string & filename) const {

  // *** WARNING:  We assume that theta_bin has size precision/32 
  // ***           (which could also be coded as precision >> 5).

  // Notation:
  // ---------
  // precision of the theta function ==> ... + O(x^(precision + 1))



  // Make the new filename: 
  // ----------------------
  string newfilename;
  if (filename == "") 
    // No filename ==> Use "THETA_DIR/boolean_theta__'QQ'__'precision'.bindata"
    newfilename = string(THETA_DIR) + "boolean_theta__" + QQ.QF_String() + "__" + MakeString(precision()) + ".bindata";
  else
    newfilename = GetAbsolutePath(filename);  



  /*
  // Diagnostic:
  cout << "boolean_theta::write -- Using filename " << newfilename << endl << endl;
  */

    
  // -------------------------------------------------------------------


  // Open the file for writing and check it opened correctly
  ofstream fileout;
  fileout.open(newfilename.c_str(),  ofstream::binary);

  // Check to make sure the file opened properly
  if (! fileout.is_open())
    { cout << "Error opening output file " << newfilename; exit(1); }

  // Write the array
  unsigned long *buffer = new unsigned long[_length()];
  copy(_theta.begin(), _theta.end(), buffer);
  fileout.write( (char*) buffer, sizeof(long) * _length());
  delete [] buffer;
  buffer = 0;

  // Close the file
  fileout.close();

}





///////////////////////////////////////////////////////////////////////
// This reads in a (Magma formatted) q-series and translates it into //
// a bitwise boolean theta function formatted as an array of long.   //
// [ Series is assumed to be increasing exponents, formatted like ]  //
// [           2 + 5*q + x^2 + ... + O(q^100)                     ]  //
// Note: We don't allow the coefficient zero.
///////////////////////////////////////////////////////////////////////

void boolean_theta::read_magma(const string & seriesfilename) {

  // Make the input filename
  string inputfilename;
  inputfilename = string(THETA_DIR) + seriesfilename; 


  /*
  // DIAGNOSTIC
  cout << inputfilename << endl;
  */


  // Open the seriesfile for reading
  ifstream seriesfile;
  seriesfile.open(inputfilename.c_str(), ios::in);
  if (! seriesfile.is_open())
    { cout << "read_magma Error: Error opening file"; exit (1); }

  
  // Initialization
  char c;
  long num;
  long pow = -1;
  bool DoneFlag = false;
    
  // Check if we're still in the current long
  do {
    
    // Check to see if we're done "O(q^...)"
    seriesfile >> c;
    if (c == 'O') {
      DoneFlag = true;
    }
    
    // Otherwise, try to read the exponent
    else {
      seriesfile.putback(c);
      
      // 1. Check to see of there is no leading coefficient (which means the coefficient is 1)
      if (c == 'q') {
	num = 1;
	c = '*';
      }
      else {
	// Get the coefficient
	seriesfile >> num;         	
	
	// Eat the "+" or "*"
	seriesfile >> c;
      }
      
      // If there's a '+', then we have a constant term
      if (c == '+') 
	pow = 0;
      else {
	
	// 2. Check if there is a linear term
	// Eat the "q +" or "q^"
	seriesfile >> c;
	seriesfile >> c;
	
	// If there's a '+', then we have a linear term
	if (c == '+') 
	  pow = 1;
	else {
	  
	  // 3. Otherwise we have a higher-order term
	  // See where to put the term if it has an exponent
	  seriesfile >> pow;
	  
	  // Eat the final "+"
	  seriesfile >> c;
	}
      }      
      
      
      // Set the bit for that power    
      if ((pow <= precision()) && (num != 0))
	set_value(pow);


      // DIAGNOSTIC
      if (pow == 10) {
	cout << "** pow = " << pow 
	     << "   c = " << c 
	     << "   num = " << num << endl;
      }
      
      
    }  // End of the  "if (c != 'O')"  clause
    
  } while ((DoneFlag == false) && (pow < precision()));
    


  // Check if the Magma series precison is smaller than the desired precision.
  if (DoneFlag == true) {
    
    // Find the current precision
    seriesfile >> c;    // Eat the '('
    //  cout << " just ate: " << c << endl;
    seriesfile >> c;    // Eat the 'q'
    //  cout << " just ate: " << c << endl;
    seriesfile >> c;    // Eat the '^'
    //  cout << " just ate: " << c << endl;
    seriesfile >> num;    // Find the precision (num - 1)
    num = num - 1; 

    
    // Error message if the magma series runs out too soon!
    if (num < precision()) {
    cout << " Error in Magma_Read: The Magma file series precision is lower than the desired precision!" <<endl;    
    cout << " The Magma file precision seems to be " << num << endl;
    cout << " while the desired precision is " << precision() << endl;
    }
    
  }
  
  
  // Close the file
  seriesfile.close();
    
}








//////////////////////////////////////////////////////////////////////
// Detailed comparison of the current boolean_ternary_theta with b. //
// Outputs the forms, levels, and first few coefficients, as well   //
// as the details around and failure of equality.                   //
//////////////////////////////////////////////////////////////////////

void boolean_theta::extended_comparison(const boolean_theta & b) {

  // Check if the forms are the same, and print the forms
  if (QQ == b.QQ)
    cout << " The quadratic forms agree, and they are: " 
	 << endl << QQ << endl;
  else
       cout << " The quadratic forms don't agree..." 
	    << endl << "  1st QF:   " << endl << QQ
	    << endl << "  2nd QF:   " << endl << b.QQ
	    << endl;
  

  // Check if the precisions are the same, and print them
  if (precision() == b.precision())
    cout << " The precisions are the same:  precision = " 
	 << precision() << endl;
  else 
    cout << " The precisions don't agree..." << endl  
	 << "    1st precision = " << precision() << endl
	 << "    2nd precision = " << b.precision() << endl;


  if ((QQ == b.QQ) && (precision() == b.precision())) {
    
    // Check if the theta series are the same.
    mpz_class i=0;
    while((i <= precision()) && (get_value(i) == b.get_value(i)))
      i++;

    // Report the results
    if (i > precision())
      cout << "The theta functions match!" << endl;
    else {
      cout << "The theta function don't match starting at i=" << i << endl;
      for (mpz_class j = max(mpz_class(i-10), 0); j <= min(mpz_class(i+20), precision()); j++)
	cout << " i = " << j << "   " 
	     << get_value(j) << "   " 
	     << b.get_value(j) << endl;      
    }
    
  }
  
}



