
// This should test the routines in count_local.cc 
// (as well as some added to main.cc)

// Warning: CURRENTLY THE COUNTLOCALTYPE ROUTINES ARE NOT TESTED...


void CountLocalTest() {

  cout << "\n\n Entering CountLocalTest: " << endl;
  cout << " ------------------------ " << endl << endl;

  Matrix_mpz Q420(4,4);

  Q420(1,1) = 1;
  Q420(2,2) = 3;
  Q420(3,3) = 5;
  Q420(4,4) = 7;

  cout << " Q420 = \n" << Q420;

  

  // Check the Matrix_mpz.EvaluateQuadratic(v,R) routine
  valarray<mpz_class> v1, v2, v3, v4;

  v1.resize(4);
  v1[0] = 1;
  v1[1] = 0;
  v1[2] = 0;
  v1[3] = 0;

  cout << "\n v1 = ";
  PrintV(v1);
  cout << endl;


  v2.resize(4);
  v2[0] = 0;
  v2[1] = 1;
  v2[2] = 0;
  v2[3] = 0;

  cout << "\n v2 = ";
  PrintV(v2);
  cout << endl;


  v3.resize(4);
  v3[0] = 0;
  v3[1] = 0;
  v3[2] = 1;
  v3[3] = 0;

  cout << "\n v3 = ";
  PrintV(v3);
  cout << endl;


  v4.resize(4);
  v4[0] = 0;
  v4[1] = 0;
  v4[2] = 0;
  v4[3] = 1;

  cout << "\n v4 = ";
  PrintV(v4);
  cout << endl;


  cout << "\n Evaluating Q420[v1] mod 5: " << Q420.EvaluateQuadratic(v1,5) << endl;
  cout << "\n Evaluating Q420[v2] mod 5: " << Q420.EvaluateQuadratic(v2,5) << endl;
  cout << "\n Evaluating Q420[v3] mod 5: " << Q420.EvaluateQuadratic(v3,5) << endl;
  cout << "\n Evaluating Q420[v4] mod 5: " << Q420.EvaluateQuadratic(v4,5) << endl;



  // Check the mpz_class ^ size_t routine
  mpz_class base;
  base = 11;
  size_t pow;

  pow = 1;
  cout << "\n\n Using mpz_class ^ size_t to evaluate: 11^1 = " << (base^pow) << endl;

  pow = 2;
  cout << "\n Using mpz_class ^ size_t to evaluate: 11^2 = " << (base^pow) << endl;

  pow = 3;
  cout << "\n Using mpz_class ^ size_t to evaluate: 11^3 = " << (base^pow) << endl;



  // Check the Increment routine
  valarray<mpz_class> vv;
  vv.resize(3);
  cout << "\n\n Looping though vectors mod 3 using Increment(vv,3): " << endl;

  for (size_t i=1; i<30; i++) {
    PrintV(vv);
    Increment(vv, 3);
  }    
  

  // Check the CountLocalNaive routine
  mpz_class total, temp;
  
  cout << "\n Using Q = Q420: " << endl;
  for(size_t i=0; i<11; i++) {
    temp = CountLocalNaive(Q420,i,11);
    total += temp;
    cout << " Number of representations of " << i << " is " << temp << endl; 
  }
  
  cout << "Total number of solutions is " << total << " and should be 11^4 = 14641 \n" << endl;
  


  // Check the IsLocalSolutionType -- still need to distinguish types 4 and 5...
  //  TO DO: This is untested for now... I'll look at it when we deal with the local densities...



  // Check the CountLocalTypeNaive routine -- This is the main routine, and all of the others are front ends...




  // Check the CountLocalType routine


  // Check the CountLocalGoodType routine


  // Check the CountLocalZeroType routine


  // Check the CountLocalBadType routine


  // Check the CountLocalBadIType routine


  // Check the CountLocalBadIIType routine

}
