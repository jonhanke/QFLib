
// Routines to read in the Eisenstein series, and compute 
// a naive lower bound using only its coefficients.




// This is a fairly naive routine to read in a q-series 
// (from Magma) and store it as a valarray<mpq_class>.

//valarray <mpq_class> ReadSeries(const string & seriesfilename) {
valarray <mpq_class> ReadSeries1(char* seriesfilename, unsigned long precision) {

  valarray <mpq_class> series;
  unsigned long SERIES_SIZE = precision;

  series.resize(SERIES_SIZE);

  ifstream seriesfile;

  seriesfile.open(seriesfilename, ios::in);

  if (! seriesfile.is_open())
    { cout << "Error opening file"; exit (1); }


  char c;
  unsigned long num;

  // Read the initial entry (q^0)
  seriesfile >> series[0];

  // Eat the "+"
  seriesfile >> c;
    //cout << " Read in " << c << endl;


  // Read the next entry (q^1)
  seriesfile >> series[1];

  // Eat the "*q +"
  seriesfile >> c;
    //cout << " Read in " << c << endl;
  seriesfile >> c;
    //cout << " Read in " << c << endl;
  seriesfile >> c;
    //cout << " Read in " << c << endl;


  for (size_t i=2; i < SERIES_SIZE; i++) {

    // Read the next entry (q^2)
    seriesfile >> series[i];
    
    // Eat the "*q^"
    seriesfile >> c;
      //cout << " Read in " << c << endl;
    seriesfile >> c;
      //cout << " Read in " << c << endl;
    seriesfile >> c;
      //cout << " Read in " << c << endl;
    
    // Check the power agrees with the index
    seriesfile >> num;
    if (i != num) 
      cout << "\n Error reading the Eisenstein series: Putting the q^" << num << " term in the " << i << "th entry." << endl;
    
    // Eat the "+"
    seriesfile >> c;
      //cout << " Read in " << c << endl;
    
  }
  
  return series;
}



////////////////////////////////////////////////////////////
// This is a more refined routine to read in a q-series   //
// (from Magma) and store it as a valarray<mpq_class>. =) //
////////////////////////////////////////////////////////////

valarray <mpq_class> ReadSeries2_half(const char* seriesfilename, long precision) {

  valarray <mpq_class> series;
  long SERIES_SIZE = precision;

  series.resize(SERIES_SIZE);

  ifstream seriesfile;

  seriesfile.open(seriesfilename, ios::in);

  if (! seriesfile.is_open())
    { cout << "ReadSeries Error: Error opening file"; exit (1); }

  char c;
  mpq_class num;
  long pow;
  bool DoneFlag = false;

  pow = -1;

  while ((DoneFlag == false) && (pow < 2*SERIES_SIZE)) {

      // cout << "\n Just finished pow = " << pow << endl;

    // Check to see if we're done "O(q^...)"
    seriesfile >> c;
      // cout << " Read in char " << c << endl;
    if (c != 'O') {
      seriesfile.putback(c);
        // cout << " Put back char " << c << endl;

      // Check to see of there is no leading coefficient (which means the coefficient is 1)
      if (c == 'q') {
	num = 1;
	c = '*';
      }
      else {
	seriesfile >> num;         	// Get the constant term
  	  // cout << " Read in num = " << num << endl;
	
	// Eat the "+" or "*"
	seriesfile >> c;
          // cout << " Read in char " << c << endl;
      }


      if (c == '+') 
	series[0] = num;
      else {
	
	if (c != '*') {
	  cout << " ReadSeries Error1: Unexpected input character " << c << " found." << endl; 
	  exit (1); 
	}
      

	  // cout << " Checking for a linear term " << endl;
	
      
	// Check if there is a linear term
	
	// Eat the "q +" or "q^"
	seriesfile >> c;
	  // cout << " Read in char " << c << endl;
	seriesfile >> c;
	  // cout << " Read in char" << c << endl;
	
	if (c == '+') 
	  series[1] = num;
	else {
	  
	  if (c != '^') {
	    cout << " ReadSeries Error2: Unexpected input character " << c << " found." << endl; 
	    exit (1); 
	  }
	  
	  
	    // cout << " Checking for higher terms " << endl;
	  
	  
	  // See where to put the term if it has an exponent
	  seriesfile >> pow;
	    // cout << " Read in pow = " << pow << endl;

	  
	  // Eat the final "+"
	  seriesfile >> c;
	    // cout << " Read in char " << c << endl;

	  
	  if ((pow > 0) && (pow < 2*SERIES_SIZE))
	    // Modification to take half of the exponent...
	    if ((pow % 2) == 1)
	      cerr << "Error in ReadSeries2_half:  The exponent " << pow << "is not even... =( " << endl;
	    else 
	      series[pow / 2] = num;
	  else {
	    /* 
	    cout << " ReadSeries Warning: Exponent " << pow << " of series is out of range (0 ... " << (SERIES_SIZE - 1) << ")." << endl; 
	    cout << " Disregarding exponent " << pow << "." << endl;
	    */
	  }
	}
      }      

    }
    else
      DoneFlag = true;
    
  }
  
  
  return series;
}



///////////////////////////////////////////
// Front end for the ReadSeries2 routine //
///////////////////////////////////////////

valarray <mpq_class> ReadSeries(const char* seriesfilename, long precision) {
  return ReadSeries2_half(seriesfilename, precision);
}







  /*
  char buffer[256];

  while (! seriesfile.eof() )
    {
      seriesfile.getline (buffer,100);
      cout << buffer << endl;
    }
  */
 



/*
template<class T>
void ReadList(istream & in, valarray<T> & v)
{
  // [ 1, 2, 3, 4 ]
  unsigned int count = 0;
  char c;
  T thingee;
  list<T> elements;
 
  // Eat the '['
  in >> c;
 
  while(true)
    {
      // Look for the ']'
      in >> c;
 
      // If it's a closing bracket, we're done looping.
      if(c == ']')
        break;
 
      // If it's a comma, we just eat it and go on.
      if(c == ',')
        continue;
                                                                                                                                                             
      // Else put it back and read a T
      in.putback(c);
      in >> thingee;
      count++;
      // Put it on the end of the list
      elements.push_back(thingee);
    }
                                                                                                                                                             
  // Having read all our thingees, we put them into a valarray
  v.resize(count);
  typename list<T>::iterator i;
  unsigned int j;
                                                                                                                                                             
  for(i = elements.begin(), j = 0; i != elements.end(); i++, j++)
    {
      v[j] = *i;
    }
}
*/
