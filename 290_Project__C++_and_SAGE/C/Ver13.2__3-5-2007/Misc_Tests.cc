



void LCM_Test() {

  // Checking the LCM routine
  cout << "The LCM of 4 and 6 is " << LCM(4,6) << endl;
  cout << "The LCM of -12 and -15 is " << LCM(-12,-15) << endl;
  cout << "The LCM of 3 and -5 is " << LCM(3,-5) << endl;

}


void Ternary_Level_Test() {
  
  // Read in the 34 ternary forms
  long QQ[34][6];
  Read34Ternaries(QQ);
  
  // Find the level of each of the 34 ternary forms
  for(long i=0; i<34; i++)
    cout << " The level of form #" << i+1 << " is " << QF_Ternary_Level(QQ[i]) << endl;

}





////////////////////////////////////////////////////////////
// Tests the ThetaPrecision routine for the 420 form and  //
// estimates how many terms we need to check the big form //
////////////////////////////////////////////////////////////

void ThetaPrecisionTest() {
  
  double Theta_precision;
  Theta_precision = ThetaPrecision(177.03, 420, 105, 7);
  cout << " The predicted Q420 theta series precision is: " << Theta_precision << endl << endl; 
  
  Theta_precision = ThetaPrecision(8232.13, 7488, 26, 1);
  cout << " The predicted theta series precision for the BIG form is: " << Theta_precision << endl << endl; 
  
}




//////////////////////////////////////
// Test code to write a binary file //
//////////////////////////////////////

void Testing_Binary_File_Operations() {

  // Make the arrays
  unsigned long x[3], y[3];
  x[0] = 12;
  x[1] = 555;
  x[2] = 99;

  // Make the filename
  char testfilename[60];  
  sprintf(testfilename, "test_long_len=%d.bin", 3); 

  // Open the file for writing and check it opened correctly
  ofstream testfileout;
  testfileout.open(testfilename, ios::out | ios::trunc | ios::binary);

  if (! testfileout.is_open())
    { cout << "Error opening output file " << testfilename; exit (1); }

  // Write the array
  testfileout.write( (char*) x, sizeof(long) * 3);

  // Close the file
  testfileout.close();

  // ---------------------------------------------------------------

  // Open the file for reading and check it opened correctly
  ifstream testfilein;
  testfilein.open(testfilename, ios::in | ios::binary);

  if (! testfilein.is_open())
    { cout << "Error opening input file " << testfilename; exit (1); }

  // Read the array
  testfilein.read( (char*) y, sizeof(long) * 3);

  // Close the file
  testfilein.close();

  // ----------------------------------------------------------------

  // Print the results
  cout << "\n\n Here are the results of our binary file read/write exercise: " << endl;
  for (long i=0; i<3; i++)
    cout << " x[" << i << "] = " << x[i] 
	 << "    y[" << i << "] = " << y[i] << endl;
  cout << endl;

}




///////////////////////////////////////////////////////
// Testing the Read/WriteTernaryThetaBinary Routines //
///////////////////////////////////////////////////////

void Testing_Binary_Theta_ReadWrite() {

  // Set the parameters
  long Q[6] = {1, 0, 0, 3, 0, 5};
  unsigned long Precision = 1000;

  // Make a theta function to use (for x^2 + 3y^2 + 5z^2)
  unsigned long Theta_vec[(Precision >> 5) + 1];
  FastShort123Binary(Theta_vec, Precision);

  // Write it
  WriteTernaryThetaBinary(Theta_vec, Q, Precision);

  // Make a new home
  unsigned long New_vec[(Precision >> 5) + 1];
  
  // Read it
  ReadTernaryThetaBinary(New_vec, Q, Precision);
  
  // Check it's the same
  for (unsigned long i=0; i < (Precision >> 5) + 1; i++)
    if (Theta_vec[i] != New_vec[i])
      cout << " Error in Testing_Binary_Theta_ReadWrite! =(    (at entry " << i << ")" << endl
	   << "    Wrote " << Theta_vec[i] << "  and read " << New_vec[i] << endl;
 
  cout << " Ok, finished checking the binary theta read/write routines now. =)" << endl;

}





////////////////////////////////////////////////////
// Tests that the 1-3-5 form encodes properly. =) //
////////////////////////////////////////////////////

void Testing_the_Encoding_Theta_Routine() {

  // Set the parameters
  long Q[6] = {1, 0, 0, 3, 0, 5};
  extern char THETA_DIR[];
  /*
  char theta_file[] = "Ternary_1-3-5__theta_series__1000.txt";
  unsigned long Precision = 1000;
  */
  char theta_file[] = "Ternary_1-3-5__theta_series__1300000.txt";
  unsigned long Precision = 1300000;


  // Encode the old (text) theta function as a new boolean (binary) theta function
  EncodeTernaryThetaBinary(theta_file, Q);


  // Read in the (text) ternary theta series
  vector<bool> text_theta;
  cout << "\n\n Started reading the ternary theta series." << endl;
  char text_filename[200]; 
  sprintf(text_filename, "%s%s", THETA_DIR, theta_file);   
  ReadSeries_Boolean(text_theta, text_filename);
  cout << " Finished reading the text ternary theta series." << endl;

  // Read in the (binary) ternary theta series
  unsigned long bin_theta[(Precision >> 5) + 1];
  cout << "\n\n Started reading the ternary theta series." << endl;
  ReadTernaryThetaBinary(bin_theta, Q, Precision);
  cout << " Finished reading the binary ternary theta series." << endl;


  // Compare them
  cout << "Comparing the two theta series: " << endl;
  for(unsigned long j=0; j < Precision; j++) {
    bool bin_bit = ((bin_theta[j >> 5] & (unsigned long) (1 << (j % 32)))) >> (j % 32);
    if (text_theta[j] != bin_bit) {
      cout << "Error in bit " << j << ". " << endl;
      cout << "  text = " << text_theta[j] << "   bin = " << bin_bit << endl;
    }
  }
  cout << "Finished comparing the two theta series. =) " << endl;

}




//////////////////////////////////////////////////////////////////////////////////
// Compares the FastShort123Binary and FastBinaryTheta routines.                //
// This tests that the FastBinaryTheta routine works, and compares their speed. //
//////////////////////////////////////////////////////////////////////////////////

void CompareThetasFor135Form() {
  
  // Set some precision constants
  const unsigned long Theta_Precision = 1000000;
  const unsigned long Theta_digits = (Theta_Precision >> 5) + 1;
  cout << "Using Precision " << Theta_Precision << " to compare the theta series," << endl;
  cout << "  which allocates " << Theta_digits << " unsigned longs." << endl;
  
  // Define the 1-3-5 form
  double Q135_double[6] = {1, 0, 0, 3, 0, 5};
  
  // Define the theta series arrays
  unsigned long theta_long[Theta_digits];
  unsigned long theta_double[Theta_digits];

  // Compute the two theta functions
  PrintTime();
  cout << " Starting FastShort123Binary " << endl;
  FastShort123Binary(theta_long, Theta_Precision);
  cout << " Finishing FastShort123Binary " << endl;
  PrintTime();
  cout << " Starting FastBinaryTheta " << endl;
  FastBinaryTheta(theta_double, Q135_double, (double) Theta_Precision);
  cout << " Finishing FastBinaryTheta " << endl;
  PrintTime();


  /*
  // Output the first 100 bits
  for(unsigned long j=0; j < 100; j++) {
    bool long_bit = ((theta_long[j >> 5] & (unsigned long) (1 << (j % 32)))) >> (j % 32);
    bool double_bit = ((theta_double[j >> 5] & (unsigned long) (1 << (j % 32)))) >> (j % 32);
    cout << "bit: " << j << "   long_bit = " << long_bit << "   double_bit = " << double_bit << "  " << endl; 
  }
  */
  

  // Compare the two theta functions
  cout << "Comparing the two theta series: " << endl;
  long ErrorCount = 0;
  for(unsigned long j=0; j < Theta_Precision; j++) {
    bool long_bit = ((theta_long[j >> 5] & (unsigned long) (1 << (j % 32)))) >> (j % 32);
    bool double_bit = ((theta_double[j >> 5] & (unsigned long) (1 << (j % 32)))) >> (j % 32);
    if ((long_bit != double_bit) && (ErrorCount < 10)){
      cout << "Error in bit " << j << ". " << endl;
      cout << "  long_bit = " << long_bit << "   double_bit = " << double_bit << " XX " << endl;
      ErrorCount++;
    }
  }
  cout << "Finished comparing the two theta series. =) " << endl;
  
}







/////////////////////////////////////////////////////////////////////////////////
// Distributed Test Code:
//
// 1. Check that MakeTernaryThetaDistributed__by_x_n agrees with 
//    MakeTernaryThetaDistributed when the slice_min is big, or the bound is small.
//




bool TestDistribuedThetaComputation() {

  /*
  // Make the "Little Methuzalah" form
  // ans set the precision
  long QQ[6] = {2, 1, -1, 4, 3, 31};
  unsigned long precision = 100000;                    // *** Note: Ideally we'd be able to pass the precision, as well as the slice size!
  unsigned long theta_length = (precision >> 5) + 1;

  // Make two vectors to store the theta function
  unsigned long theta1[theta_length];
  unsigned long theta2[theta_length];


  // Make and read in the two forms
  MakeTernaryTheta(QQ,precision);
  ReadTernaryThetaBinary(theta1, QQ, precision);

  MakeTernaryThetaDistributed__by_x_n(QQ,precision);
  ReadTernaryThetaBinary(theta2, QQ, precision);

  
  // Compare them
  bool equal_flag = true;
  for(unsigned long i=0; i < theta_length; i++)
    if (theta1[i] != theta2[i]) {
      equal_flag = false;
      cout << "Error in TestDistribuedThetaComputation:" << endl;
      cout << "  The two little methuzalah theta functions differ in byte " << i << endl;
      cout << "    Undistributed: " << theta1[i] << "   Distributed: " << theta2[i] << endl;
      cout << "  Here are the previous entries: " << endl;
      cout << "    Undistributed: " << theta1[i-1] << "   Distributed: " << theta2[i-1] << endl;


    }

  return equal_flag;
  */      

  return false;


}




