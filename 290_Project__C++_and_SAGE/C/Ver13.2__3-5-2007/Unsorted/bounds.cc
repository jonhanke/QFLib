

///////////////////////////////////////////////////////////////////////
// Compute the local constants for each form and save them to a file //
///////////////////////////////////////////////////////////////////////

void FindLowerBounds() {

  int ListEnd = 6560;

  // Select a range of forms to compute
  cout << " Which form (>=1 and <= " << ListEnd << ") would you like to start with? ";
  size_t begin_form;
  cin >> begin_form;

  cout << " Which form (>=" << begin_form << " and <= " << ListEnd << ") would you like to end with? ";
  size_t end_form;
  cin >> end_form;

  assert(1 <= begin_form <= end_form <= ListEnd);


  valarray <Matrix_mpz> form_list;
  valarray <mpz_class> lvl_list;
  valarray <float> cusp_const_list;

  valarray <mpq_class> local_const_list;
  valarray <float> computed_lower_bound_list;
  valarray <mpq_class> approximate_lower_bound_list;

  form_list.resize(ListEnd);
  lvl_list.resize(ListEnd);
  cusp_const_list.resize(ListEnd);

  local_const_list.resize(ListEnd);
  computed_lower_bound_list.resize(ListEnd);
  approximate_lower_bound_list.resize(ListEnd);


  // Read in the 6560 quadratic forms 

  cout << "\n Reading the 290 cusp data file " << endl;
  Read290Data("/home/postdoc/jonhanke/290_Project/290_Cusp_info/290-cusp-all.m", form_list, lvl_list, cusp_const_list);
  cout << "\n Finished reading the 290 cusp data file \n" << endl;





  // Loop through several forms to check their local factors...  
 
  for(size_t k = begin_form-1; k < end_form; k++) {
    

    // Redirect cout to an information file for each form
    streambuf *psbuf;
    ofstream filestr;
    stringstream Form_output_filename_stream;
    Form_output_filename_stream << "/home/postdoc/jonhanke/290_Project/F4_Complete_Output/Form" << k+1 << "_F4_output.txt";
    const char* Form_output_filename;
    Form_output_filename = Form_output_filename_stream.str().c_str();
    filestr.open (Form_output_filename);

    psbuf = filestr.rdbuf();
    cout.rdbuf(psbuf);


    // Start processing the (k+1)st form
    Matrix_mpz Q;
    valarray <mpq_class> Eis;

    Q = form_list[k];
    
    cout << "\n\n Looking at form #" << k+1 << ":\n" << endl;
    //    cout << Q << "\n\n";


    // Find the Ansotropic primes
    valarray <mpz_class> SS;
    SS = AnisotropicPrimes(Q);
    cout << " The anisotropic primes (for form #" << k+1 << ") are: " << SS << endl;
    
    
    // Prepare to read in the Eisenstein seies
    stringstream filename; 
    const char* Eisfile;
    filename << "/home/postdoc/jonhanke/290_Project/EIS_DATA/tmp/Jon_Eis/Form" << k+1 << "_Eis_10000.txt";
    Eisfile = filename.str().c_str();


    // The precision can be up to 10,000 for this Eisenstein series
    Eis = ReadSeries(Eisfile, 10000);
    
    cout << " Read in the Eisenstein series vector " << endl;
      // cout << Eis << endl << endl;

    

    //    /*

    // Prepare to write the F4_summary_file
    stringstream F4_summary_filename; 
    const char* F4summaryfile;
    F4_summary_filename << "/home/postdoc/jonhanke/290_Project/F4_Summaries/Form" << k+1 << "_F4_summary.txt";
    F4summaryfile = F4_summary_filename.str().c_str();


    ofstream outsfile;
    outsfile.open(F4summaryfile, ios::out);

    if (! outsfile.is_open())
      { cout << "Error opening file"; exit (1); }


    outsfile << "\n\n Looking at form #" << k+1 << ":\n" << endl;



    // Compute the theoretical local constant Lambda_4^{\hat}
    mpq_class local_const;
    local_const = new_C4_rational(Q);

    // FIX THIS!!!
    // --> Quick fix: Extra factor of 4 (only for n=4) since we altered the local normal forms from 2Q --> Q
    local_const *= 4;


    // Compute the numerical (approximate) local constant Lambda_4^{\hat} from the first few Eisenstein coefficients
    mpq_class num_local_const;
    num_local_const = Lambda4_test(Q, Eis);


    // Check that the numerical local constant is  >=  the theoretical one
    if (num_local_const >= local_const) { 
      outsfile << " Ok, the numerical bound is " << num_local_const << " and the theoretical bound is " << local_const << endl;
      cout << " Ok, the numerical bound is " << num_local_const << " and the theoretical bound is " << local_const << endl;
    }
    else {
      outsfile << " ERROR! The numerical bound is " << num_local_const << " and the theoretical bound is " << local_const << endl;
      cout << " ERROR! The numerical bound is " << num_local_const << " and the theoretical bound is " << local_const << endl;
      exit(0);
    }



    // Compute the lower bound for F_4 to ensure representability by the form
    float F4_lower_bound;
    F4_lower_bound = cusp_const_list[k] / local_const.get_d();
    
    outsfile << " The cusp constant is  " << cusp_const_list[k] << endl;
    outsfile << " The local constant is " << local_const << endl;
    outsfile << " This gives the lower bound for F_4 as: " << F4_lower_bound << endl << endl << endl << endl; 

    cout << " The cusp constant is  " << cusp_const_list[k] << endl;
    cout << " The local constant is " << local_const << endl;
    cout << " This gives the lower bound for F_4 as: " << F4_lower_bound << endl << endl << endl << endl; 

    outsfile.close();

//    */


    //    /*


    //    float F4_lower_bound;
    //    F4_lower_bound = 3.14159265;

    //    cout << F4_lower_bound << endl;


    // Write all of this to datafile(s)
    ////////////////////////////////////

    // Prepare to write the F4_constant_file
    stringstream F4_filename; 
    const char* F4file;
    F4_filename << "/home/postdoc/jonhanke/290_Project/F4_Constants/Form" << k+1 << "_F4_bound.txt";
    F4file = F4_filename.str().c_str();

    //    cout << " F4file is " << F4file << endl;

    ofstream outfile;
    outfile.open(F4file, ios::out);

    if (! outfile.is_open())
      { cout << "Error opening file"; exit (1); }


    // Write the constant
    outfile << F4_lower_bound << endl;


    outfile.close();


    // Close the complete output file for this form
    filestr.close();

    //    */

  } 
  

} 


/////////////////////////////
// Find the next prime > x //
/////////////////////////////
mpz_class NextPrime(mpz_class x) {

  mpz_t y1;
  mpz_init(y1);


  // Get the next prime
  mpz_nextprime(y1, x.get_mpz_t());


  // Copy the answer into an mpz_class	     
  mpz_class y;
  y = mpz_class(y1);
  
  
  // Clean up the mpz_t variables
  mpz_clear(y1); 

  return y;
}




//////////////////////////////////////////////////////////
// Find a list of all primes allowed to have F4(p) <= X //
//////////////////////////////////////////////////////////

void AllowedPrimesList(mpq_class X, mpz_class N, mpz_class char_top) {


  // Step through all primes (smallest to largest) to find the product with the most factors.

  mpq_class X2;
  X2 = X*X;

  mpq_class F4_Prod2;
  F4_Prod2 = 1;

  mpq_class m_Prod;
  m_Prod = 1;

  mpz_class p;
  p=2;


  
  
  size_t count;
  count = 0;

  while (F4_Prod2 <= X2) {
    F4_Prod2 = F4_Prod2 * mpq_class(p,4);
    if ((p % N != 0) && (KroneckerSymbol(char_top, p) == -1))
      F4_Prod2 = F4_Prod2 * mpq_class(p-1, p+1) * mpq_class(p-1, p+1);
    m_Prod = m_Prod * p;
    p = NextPrime(p);
    count++;
  }


  cout << " The prime count is: " << count << endl;
  cout << " The largest prime is: " << p << endl;
  cout << " The product of the F_4's squared is: " << F4_Prod2 << endl;
  cout << " The product of the primes is: " << m_Prod << endl;
  cout << " X is " << X << endl;
  cout << " X^2 is " << X2 << endl;

}
