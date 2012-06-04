

/////////////////////////////////////////////////////////////////
// Routine to analyze the local densities at those numbers m   //
// where the product formula and Eisenstein series don't match //
/////////////////////////////////////////////////////////////////

void AnalyzeDensities290() {

  int EisMax = 100;
  size_t Starting_Form = 1;  // Must be >= 1
  size_t List_End = 6560;

  valarray <Matrix_mpz> form_list;
  valarray <mpz_class> lvl_list;
  valarray <float> cusp_const_list;

  form_list.resize(List_End);
  lvl_list.resize(List_End);
  cusp_const_list.resize(List_End);


  // Read in the 6560 quadratic forms 

  cout << "\n Reading the 290 cusp data file " << endl;
  Read290Data("/home/postdoc/jonhanke/290_Project/290_Cusp_info/290-cusp-all.m", form_list, lvl_list, cusp_const_list);
  cout << "\n Finished reading the 290 cusp data file \n" << endl;


  // Choose to list or analyze the discrepancies
  char ch;
  bool analyze_flag;
  
  do {
    cout << " Analyze or list the Eisenstein/product formula discrepancies? (a/l) "; 
    cin >> ch;    
  } while ((ch != 'a') && (ch != 'l'));

  if (ch == 'a') 
    analyze_flag = true;
  else
    analyze_flag = false;


  // Loop through several forms to check their local factors...
 
  Matrix_mpz Q;
  PowerSeries<mpq_class> Eis;
  mpz_class m;

  
  for(size_t k = (Starting_Form - 1) ; k < List_End; k++) {
  //  for(size_t k=6460; k<6461; k++) {
    
    Q = form_list[k];
    
    cout << "\n\n Looking at form #" << k+1 << ":\n" << endl;
    cout << Q << "\n\n";

    //    /*

        cout << " It's level is " << Q.QFLevel() << endl << endl;
    
    valarray <mpz_class> SS;
    SS = Q.AnisotropicPrimes();
        cout << " The anisotropic primes (for form #" << k+1 << ") are: " << SS << endl;


    // /*
	// Test for anisotropicity of p=31 for Form #6561
	cout << " LocalDiagonal(Q, 31) = " << Q.LocalDiagonal(31) << endl;
	cout << " Q.Determinant() = " << Q.Determinant() << endl;
	cout << " IsPadicSquare(Q.Determinant(), 31) = " << IsPadicSquare(Q.Determinant(), 31) << endl;
	cout << " HasseInvariant(Q, 31) = " << Q.HasseInvariant(31) << endl;
	cout << " HilbertSymbol(-1, -Q.Determinant(), 31) = " << HilbertSymbol(-1, -Q.Determinant(), 31) << endl;


    // */


    stringstream filename; 
    const char* Eisfile;
    filename << "/home/postdoc/jonhanke/290_Project/EIS_DATA/tmp/Jon_Eis/Form" << k+1 << "_Eis_10000.txt";
    Eisfile = filename.str().c_str();


    // The precision can be up to 10,000 for this Eisenstein series
    Eis = Q._GetEisData(Eisfile, 10000);  // Using the custom filename, but the default precision. 
    //  This is the old Eisenstein series routine
    //    Eis = ReadSeries(Eisfile, EisMax);
    
        cout << " Read in the Eisenstein series vector " << Eis << endl << endl;

        cout << " The Local Normal form of form #" << k+1 << " at 2 is: \n" << Q.GetLocalNormal(2) << endl; 

	//	cout << "  What's wrong? " << endl;    
       

	/*
    for (m=1; m<Eis.size(); m++) {
      cout << " m = " << m << "  ==>   Prod = " << SiegelProduct(Q,m) << "   Eis = " << Eis[m.get_ui()] << endl;
    }
	*/


    for (m=1; m<=Eis.Precision(); m++) {
       cout << " Checking the coefficient at m = " << m << ": ";
       if (Q.CheckSiegel(m, Eis) == true) {
	 cout << "";
	 cout << " ok " << endl; 
       }
       else {
	 cout << " trouble... " << endl;
	 if (analyze_flag == true)	  
	   break;
      }
    }



        cout << " The Local Normal form of form #" << k+1 << " at 2 is: \n" << Q.GetLocalNormal(2) << endl; 

        cout << " The Local Normal form of form #" << k+1 << " at 3 is: \n" << Q.GetLocalNormal(3) << endl; 

    //    */


    // Analyze the first discrepancy in local densities...
    if ((analyze_flag == true) && (m < EisMax)) {
      
      valarray <mpz_class> PrimeList;
      PrimeList = PrimeDivisors(Q.QFLevel());
      
      for (size_t l=0; l<PrimeList.size(); l++) {
	cout << " Analyze the local factor at p = " << PrimeList[l] << "? (y/n) ";
	cin >> ch;
	if (ch != 'y') 
	  cout << endl << endl;
	else {
	  Q.CheckLocalDensity(PrimeList[l], m);
	}
	
	cout << "\n\n ==================================== \n\n";
	
      }
      
      
      mpz_class p;
      
      ch = 'y';
      
      while (ch == 'y') {
	cout << " Would you like to analyze the local density at any other primes and numbers? (y/n) ";
	cin >> ch;
	cout << endl << endl;
	
	
	// Add support here for analyzing other primes and numbers too...
	if (ch == 'y') {
	  cout << " Which prime would you like to analyze? ";
	  cin >> p;
	  cout << endl << endl;
	  
	  cout << " Which number would you like to analyze? ";
	  cin >> m;
	  cout << endl << endl;
	  
	  Q.CheckLocalDensity(p, m);
	}
      }   
    }
    else 
      if (m < 10)
	exit(0);


  }  

}
