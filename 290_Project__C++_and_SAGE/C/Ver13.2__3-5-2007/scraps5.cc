

//////////////////////////////////////////
// Test the modified local density routines
//////////////////////////////////////////

void TestModifiedDensities() {

  vector<mpz_class> primes;
  primes.push_back(2);
  primes.push_back(3);
  primes.push_back(5);
  primes.push_back(7);
  primes.push_back(11);
  primes.push_back(13);
  primes.push_back(17);
  primes.push_back(19);
  primes.push_back(23);
  primes.push_back(29);
  primes.push_back(31);
  primes.push_back(37);
  primes.push_back(41);

  for(long i=0; i<primes.size(); i++) {
    
    // Make the local normal matrix [1, 1, 2, 58] at m=1 and p=29
    Matrix_mpz Q(4,4);
    Q(1,1) = 1;
    Q(2,2) = 1;
    Q(3,3) = 2;
    Q(4,4) = 58;
    
    // Making m and p
    mpz_class p;
    p = primes[i];
    mpz_class m;
    m = 1;
    
    // Finding the local density
    cout << " Using the Matrix: " << endl << Q;
    cout << " Now finding its local density at p = " << p << "  and m = " << m << endl;
    cout << "   Local Density = " << Q.Local_Density(p, m) << endl;
    
    
    
    // Make the local normal matrix 2*[1, 1, 2, 58] at m=1 and p=29
    Matrix_mpz Q1(4,4);
    Q1(1,1) = 2;
    Q1(2,2) = 2;
    Q1(3,3) = 4;
    Q1(4,4) = 116;
    
    // Finding the local density
    cout << " Using the Matrix: " << endl << Q1;
    cout << " Now finding its local density at p = " << p << "  and m = " << m << endl;
    cout << "   Local Density = " << Q1.Local_Density(p, m) << endl;
    
  }

}



