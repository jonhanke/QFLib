
/////////////////////////////////////////////////////////
// Get the precision of a series (in the Magma format) //
/////////////////////////////////////////////////////////

mpz_class GetMagmaSeriesPrecision(const char* seriesfilename) {  

// TO DO: Need to create an error message if the file ends before reading the 'O' character...

  ifstream seriesfile;
  seriesfile.open(seriesfilename, ios::in);
  if (! seriesfile.is_open())
    { cout << "GetMagmaSeriesPrecision Error: Error opening file"; exit (1); }

  // Read characters until we hit 'O(q^...)'
  char c;
  seriesfile >> c;
  while (c != 'O')
    seriesfile >> c;

  // Find the exponent of the error term
  seriesfile >> c;  // Read '('
  seriesfile >> c;  // Read 'q'
  seriesfile >> c;  // Read '^'

  mpz_class pow;
  seriesfile >> pow;  // Read the exponent

  return(pow);
}




//////////////////////////////////////////////////////////////
// This is a more refined routine to read in a q-series     //
// (from Magma) and store it as the vector<bool> series. =) //
//////////////////////////////////////////////////////////////

void ReadSeries_Boolean(vector<bool> & series, const char* seriesfilename) {
  
  // Get the series precision and resize the return_list accordingly
  series.resize(0);
  mpz_class Series_precision;
  Series_precision = GetMagmaSeriesPrecision(seriesfilename);
  if (!Series_precision.fits_ulong_p()) {
    cout << " Error in ReadSeries_Boolean: The series precison " << Series_precision 
	 << " doesn't fit in an unsigned long... =(" << endl;
    exit(1);  // Exit the program
  }
  series.resize(Series_precision.get_ui(), false);  // Resize the vector and make it empty.
  
  
  // Open the seriesfile for reading
  ifstream seriesfile;
  seriesfile.open(seriesfilename, ios::in);
  if (! seriesfile.is_open())
    { cout << "ReadSeries Error: Error opening file"; exit (1); }
  
  
  
  // Read in the series (as a boolean vector)
  char c;
  mpq_class num;
  long pow = -1;
  bool DoneFlag = false;
  
  while ((DoneFlag == false) && (pow < Series_precision)) {
    
    // Check to see if we're done "O(q^...)"
    seriesfile >> c;
    if (c != 'O') {
      seriesfile.putback(c);
      
      // Check to see of there is no leading coefficient (which means the coefficient is 1)
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
      if ((c == '+') && (num != 0)) 
	series[0] = true;
      else {
	
	// ERROR CHECKING: We expect a '*' here.
	if (c != '*') {
	  cout << " ReadSeries_boolean Error1: Unexpected input character " << c << " found." << endl; 
	  exit (1); 
	}
	
	
	// Check if there is a linear term
	//       cout << " Checking for a linear term " << endl;
	
	// Eat the "q +" or "q^"
	seriesfile >> c;
	seriesfile >> c;
	
	// If there's a '+', then we have a linear term
	if ((c == '+') && (num != 0)) 
	  series[1] = true;
	else {
	  
	  // Otherwise we have a higher-order term
	  if (c != '^') {
	    cout << " ReadSeries_boolean Error2: Unexpected input character " << c << " found." << endl; 
	    exit (1); 
	  }
	  
	  
	  //     cout << " Checking for higher terms " << endl;
	  
	  
	  // See where to put the term if it has an exponent
	  seriesfile >> pow;
	    // cout << " Read in pow = " << pow << endl;
	  
	  // Eat the final "+"
	  seriesfile >> c;
	    // cout << " Read in char " << c << endl;

	  
	  if ((pow > 0) && (pow < Series_precision) && (num != 0))
	    series[pow] = true;
	  else {	    
	    cout << " ReadSeries_boolean Warning: Either the exponent " << pow 
		 << " of series is out of range (0 ... " << (Series_precision - 1) << ")," 
		 << " or the coefficient " << num << " is zero." << endl; 
	    cout << " Disregarding exponent " << pow << "." << endl;
	    
	  }
	}
      }      
      
    }
    else
      DoneFlag = true;
    
  }
  
  //  return series;
}


