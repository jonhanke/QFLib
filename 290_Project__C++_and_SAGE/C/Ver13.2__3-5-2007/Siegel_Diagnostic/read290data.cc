
// Reads the datafile 290-cusp-all.m
//valarray <Matrix_mpz> Read290Data(char* filename, valarray<mpz_class> level_data, valarray <float> cusp_constants) {
void Read290Data(char* filename, valarray <Matrix_mpz> & form_data, valarray<mpz_class> & level_data, valarray <float> & cusp_data) {
//void Read290Data(char* filename) {

  /*
  valarray <mpq_class> series;
  long SERIES_SIZE = precision;

  series.resize(SERIES_SIZE);
  */

  ifstream datafile;

  datafile.open(filename, ios::in);

  if (! datafile.is_open())
    { cout << "Error opening file"; exit (1); }


  char c;
  int num;
  mpz_class lvl;
  mpq_class temp_entry;
  Matrix_mpz M(4,4);
  float cusp_const, error_const;



  // Eat the initial "data := ["
  datafile >> c;
    //cout << " Read in " << c << endl;
  datafile >> c;
    //cout << " Read in " << c << endl;
  datafile >> c;
    //cout << " Read in " << c << endl;
  datafile >> c;
    //cout << " Read in " << c << endl;
  datafile >> c;
    //cout << " Read in " << c << endl;
  datafile >> c;
    //cout << " Read in " << c << endl;
  datafile >> c;
    //cout << " Read in " << c << endl;


  // Read in all of the data... =)
  for(size_t k=0; k<6560; k++) {
  //  for(size_t k=0; k<20; k++) {

    /*
    // Skip some lines
    cout << "\n\n";
    cout << " Starting with k = " << k << endl;
    */

    // Eat the initial "<"
    datafile >> c;
      //cout << " Read in " << c << endl;
    
    // Read in the number of the quadratic form
    datafile >> num;
    //    cout << " num = " << num << endl;

    // Eat the ","
    datafile >> c;
      //cout << " Read in " << c << endl;
  

    // Read in the Matrix
    for(size_t i=1; i<=4; i++)
      for(size_t j=1; j<=4; j++) {
	
	// Eat the "[" or ","
	datafile >> c;
          //cout << " Read in " << c << endl;

	// Read in the number of the quadratic form
	datafile >> temp_entry;
	M(i,j) = 2*temp_entry;
	
      }

    form_data[num-1] = M;

    /*
    cout << "The matrix is \n" << endl;
    cout << M << endl;
    */
    
    // Eat the closing "],"
    datafile >> c;
      //cout << " Read in " << c << endl;
    datafile >> c;
      //cout << " Read in " << c << endl;


    // Read the level
    datafile >> lvl;
    level_data[num-1] = lvl;
    //    cout << " Level = " << lvl << endl;

    // Eat the ","
    datafile >> c;
      //cout << " Read in " << c << endl;


    // Read the cusp constant
    datafile >> cusp_const;
    cusp_data[num-1] = cusp_const;
    //    cout << " The cusp constant is " << cusp_const << endl;

    // Eat the ","
    datafile >> c;
      //cout << " Read in " << c << endl;
  

    // Read the error term
    datafile >> error_const;
    //    cout << " The error constant is " << error_const << endl;

    // Eat the ">," or ">]"
    datafile >> c;
      //cout << " Read in " << c << endl;
    datafile >> c;
      //cout << " Read in " << c << endl;
    
  }


  // Eat the final ";"
  datafile >> c;
    //cout << " Read in " << c << endl;

  
}
