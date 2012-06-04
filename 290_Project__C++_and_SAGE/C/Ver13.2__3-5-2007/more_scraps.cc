
////////////////////////////////////////////////////
// Tests some simple array and bitwise operations //
////////////////////////////////////////////////////

void  MiscellaneousArrayAndBitOperationTests() {
  cout << "\n sizeof(long) = " << sizeof(long) << endl;


  unsigned long test[5];
  for (long i=0; i<5; i++)
    test[i] = 0;


  cout << "\n\n This is test[5]" << endl;
  for (long i=0; i<5; i++)
    cout << test[i] << " ";

  test[0] = test[0] ^ 2;  

  cout << "\n\n This is test[5]" << endl;
  for (long i=0; i<5; i++)
    cout << test[i] << " ";

  test[0] = test[0] ^ 6;  

  cout << "\n\n This is test[5]" << endl;
  for (long i=0; i<5; i++)
    cout << test[i] << " ";

  cout << endl << endl;
}



////////////////////////////////////////////////////////
// Makes a list of universal primes...                //
//   (shows there is too much error for small primes, //
//   which is why we made the optimized routine.)     //
////////////////////////////////////////////////////////

void MakingUniversalPrimesListTesting() {

  vector<long> Plist;
  vector<double> F4list;

  vector<mpz_class> Tlist;


  // Make the universal list of primes...
  cout << "\n\nStarting to list the universal eligible primes at: " ;
  PrintTime();
  PrimeListFromBoundUniversal(Plist, F4list, B, all_primes);
  cout << " There are " << Plist.size() << " eligible prime numbers." << endl;
  cout << "Finished listing the universal eligible primes at: " ;
  PrintTime();

  cout << " The last 7 eligible primes and their F4(p) values are: " << endl;
  PrintHeadV(Plist, 7);
  PrintTailV(Plist, 7);
  PrintHeadV(F4list, 7);
  PrintTailV(F4list, 7);
}



////////////////////////////////////////////////////////////////////////
// Code to make lists of eligible squarefree numbers in various ways. //
//   (All of these are superceded by the incremental routine          //
//   which makes and checks numbers at the same time.                 //
////////////////////////////////////////////////////////////////////////

void MakingEligibleSquareFreeLists() {  

  // Make a list of eligible square-free numbers
  /*
  FasterSharpList(Tlist, B, Plist, F4list);
  cout << "\n There are " << Tlist.size() << " eligible square-free numbers." << endl;
  PrintTime();
  */


  // Make a list of eligible square-free numbers
  /*
  SharpList(Tlist, B, Plist, F4list);
  cout << "\n There are " << Tlist.size() << " eligible square-free numbers." << endl;
  PrintTime(); 
  */


  // Make a list of eligible square-free numbers
  /*
  ForgetfulSharpList(Tlist, B, Plist, F4list);
  cout << "\n There are " << Tlist.size() << " eligible square-free numbers." << endl;
  PrintTime();
  */


  // Make a list of eligible square-free numbers
  /*
  LazySharpList(Tlist, B, Plist, F4list);
  cout << "\n There are " << Tlist.size() << " eligible square-free numbers." << endl;
  PrintTime(); 
  */

}




///////////////////////////////////////////////////
// Tests the speed of the various routines which //
// make the x^2 + 3y^2 + 5z^2 theta function.    //
///////////////////////////////////////////////////

void MakingTernary135ThetaTesting() {
  
  vector<mpz_class> new_theta;

  // Short123 Running Times:
  // =======================
  // - 100,000 ==> 29 seconds
  // - With x^(3/2) scaling, 1,000,000 should take about 31.6 times as long, which is about 15 minutes.
  // - This means that the 10,000,000 with take about 8 hours!?!... =(

  /*
  PrintTime();
  Short123(new_theta, mpz_class(100000));
  PrintTime();

  cout << endl << endl;
  */

  // FastShort123 Running Times:
  // ===========================
  // - 100,000 ==> 4 seconds  =)     <== This is much faster! =)


  PrintTime();
  FastShort123(new_theta, 100000);
  PrintTime();


  /*
  // Write the theta function
  cout << "Here's my 1-3-5 theta function: " << endl;
  for(unsigned long j=0; j < new_theta.size(); j++) 
    cout << j << ":  " << new_theta[j] << endl;
  */




  // NOTE: We must declare the size of the vector so it doesn't get erased when we leave the function!
  unsigned long newV1[(1000000 >> 5) + 1];
  PrintTime();
  FastShort123Binary(newV1, 1000000);
  PrintTime();

  unsigned long newV2[(10000000 >> 5) + 1];
  PrintTime();
  FastShort123Binary(newV2, 10000000);
  PrintTime();


  //  /*
  // Output the binary theta function... =)
  cout << "Here's my binary 1-3-5 theta function: " << endl;
  for(unsigned long j=0; j < 2; j++) {
    cout << "The long is: " << newV1[j] << endl;
    for(unsigned long k=0; k < 32; k++) 
      cout << (32*j + k) << ":  " << ((newV1[j] & (1 << k)) == (unsigned long) (1 << k)) 
	   << "    " << k << " mod 32 = " << (k % 32) << endl;
  }
  //  */
  
}  




////////////////////////////////////////////////////////
// Checks the squarefree list for the 420 form by     //
// using the textfile for the ternary theta function. //
////////////////////////////////////////////////////////

void Check420SquarefreeByTextTheta() {

  // Read in the text ternary theta series
  vector<bool> Ternary_series;

  cout << "\n\n Started reading the ternary theta series." << endl;
  PrintTime();
  char text_filename[200]; 
  sprintf(text_filename, "%s%s", THETA_DIR, theta_file);   
  ReadSeries_Boolean(Ternary_series, text_filename);
  cout << " Finished reading the ternary theta series." << endl;
  PrintTime();

  cout << "Here are the first 50 non-zero terms of the ternary theta series: " << endl;
  for(unsigned long j=0; j < 50; j++) 
    if (Ternary_series[j] == true)
      cout << j << ", ";
  cout << endl << endl;    


  // Check the squarefree numbers
  NewSharpList(Tlist, B, Plist, F4list, Ternary_series, 7);
}





////////////////////////////////////////////////////////
// Tests the old CheckNumbers routine on the 420 form //
////////////////////////////////////////////////////////

void TestingCheckNumbersRoutineOn420() {

  // Read in the ternary theta series
  vector<bool> Ternary_series;

  cout << " Started reading the ternary theta series." << endl;
  PrintTime();
  ReadSeries_Boolean(Ternary_series, theta_file);
  cout << " Finished reading the ternary theta series." << endl;
  PrintTime();

  cout << "Here are the first 50 non-zero terms of the ternary theta series: " << endl;
  for(unsigned long j=0; j < 50; j++) 
    if (Ternary_series[j] == true)
      cout << j << ", ";
  cout << endl << endl;    


  // Check them against the ternary theta series
  vector<mpz_class> Exception_List;
  long Checking_depth = 5;
  Exception_List = CheckNumbers(Tlist, Checking_depth, Ternary_series, 7);
  //      Exception_List = CheckNumbers(Test_Tlist, Checking_depth, Ternary_series, 7);

  cout << "There are " << Exception_List.size() << " exceptions. " << endl;

  // Write the list of exceptions...
  cout << "They are: " << endl;
  cout << "[ ";
  for(unsigned long j=0; j < (Exception_List.size()-1); j++) 
    cout << Exception_List[j] << ", ";
  cout << Exception_List[Exception_List.size()-1] << "] ";



  /*
  // TESTING: Check to see that 2 and 22 are in the Tlist;
  for(unsigned long j=0; j < Tlist.size(); j++) 
    if ((Tlist[j] == 2) || (Tlist[j] == 22))
      cout << " Found " << Tlist[j] << " at entry " << j << endl;
  cout << endl;


  // TESTING: Test code for CheckNumbers
  vector<mpz_class> Test_Tlist;
  Test_Tlist.push_back(mpz_class("6469693230"));
  Test_Tlist.push_back(150);
  Test_Tlist.push_back(2);
  Test_Tlist.push_back(22);

  vector<mpz_class> Exception_List;
  long Checking_depth = 5;
  Exception_List = CheckNumbers(Test_Tlist, Checking_depth, Ternary_series, 7);
  */
}
