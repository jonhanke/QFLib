
/*  Be sure to make a header file for this...
---------------------------------------------
forward C4_squareclass_constants_odd_rational;
forward C4_squareclass_constants_even_rational;
*/


///////////////////////////////////////////////////////////
// Make a vector of primitive square class repn's in Z_p //
///////////////////////////////////////////////////////////

valarray<mpz_class> PrimSqRepns(mpz_class p) {

  // Should check if p is prime...
  // cout << "Error in UnitSqRepns:  p is not prime!";

  valarray<mpz_class> vec;
                                                                                                                   
  if (p == 2) {

    //    return [1,3,5,7,2,6,10,14];
    vec.resize(8);
    vec[0] = 1;
    vec[1] = 3;
    vec[2] = 5;
    vec[3] = 7;
    vec[4] = 2;
    vec[5] = 6;
    vec[6] = 10;
    vec[7] = 14;
  }
  
  else {

    //    return [1, NonResidue(p), p, p*NonResidue(p)];
    vec.resize(4);
    vec[0] = 1;
    vec[1] = NonResidue(p);
    vec[2] = p;
    vec[3] = p * NonResidue(p);
  }
  
  return vec;
}



//////////////////////////////////////////////////////////////////////////////////////
// Takes in a list of representatives for the primitive squareclasses in Z_p and    //
// returns the minimum of beta_p(m) and beta_p(T') C_p(T') within each square class //
//////////////////////////////////////////////////////////////////////////////////////

valarray<mpq_class> C4_squareclass_constants(Matrix_mpz Q, mpz_class p, valarray<mpz_class> sqclasslist) {
  
  mpz_class N;
  N = Q.QFLevel();
  // cout << "Level = " << N << endl;
  
  size_t list_length;
  list_length = sqclasslist.size();
  
  valarray<mpq_class> const_list;
  const_list.resize(list_length);
  
  mpz_class t;
  unsigned long Stable;
  
  for (size_t i=0; i<list_length; i++) {
    
    t = sqclasslist[i];
    
    Stable = ((Valuation(N, p) + 1) / 2); 	// In fact, this is bigger than we need -- since it would suffice
                                                // without multipling by T below, but to make it exact is more work
                                                // with no clear benefit.  
                                                // (Note: This is forced to be an integer by integer division.)
    
    //  density_vector := [ LocalDensity(Q, p, t * p^(2*i)): i in [0..Stable+1] ];  // Should be primitively constant from Stable+1 onwards! =)
    //  prim_density_vector := [ LocalDensityGood(Q, p, t * (p^(2*i))) + LocalDensityBad(Q, p, t * (p^(2*i))) : i in [0..Stable+3] ];
    //  print "For the square-class ", t, "Z_", p, " we have the local density vector ";
    //  print density_vector;
    //  print "   we also have the local primitive density vector ";
    //  print prim_density_vector;
    //  print "";
    
    /*
      cout << "\n Starting to create the density vector for p = " << p << "\n" << endl;
      cout << " p = " << p << "  and  N = " << N << endl;
      cout << " Using the square-class entry t = " << t << endl;
      cout << " Stable = " << Stable << endl << endl;
    */
    
    // Assign the product C_p beta_p(T') to the last entry
    
    valarray <mpq_class> density_vector;
    density_vector.resize(Stable+3);

    
    Matrix_mpz Q_normal;
    Q_normal = LocalNormal(Q, p);
    
    
    if (IsIsotropic(Q, p) == true) {
      
      for (size_t j=0; j<=Stable+2; j++) 
	density_vector[j] = Local_Density(Q_normal, p, t * (p^(2*j)));  // Should be primitively constant from Stable+1 onwards! =)
      

      mpq_class Cp, Cp1;
      Cp = ((density_vector[Stable+2] / density_vector[Stable+1]) * (p*p) - 1) / ((p*p) - 1);   
      Cp1 = Minimum(Cp, 1);


      // Check that we don't have any anisotropic primes by mistake!
      assert (Cp1 > 0);

      
      // /*
      cout << " Stable = " << Stable << endl; 
      cout << " density_vector = " << density_vector << endl;
      cout << " C_p(T') = " << Cp << endl;
      cout << " C'_p(T') = " << Cp1 << endl << endl;
      // */
      
      density_vector[Stable+2] =  Cp1 * density_vector[Stable+1];
    }  
    else {
      
      // Note: density_vector should be primitively constant from Stable+1 onwards! =)
      for (size_t j=0; j<=Stable+2; j++) 
	density_vector[j] = Local_Density(Q_normal, p, t * (p^(2*j))) * (p^(Valuation(t,p) + 2*j));
      //      density_vector[Stable+2] =  density_vector[Stable+1];
      
    }
    
    
    // Find the smallest non-zero entry in density_vector
    mpq_class t_min;
    t_min = 0;
    for (size_t j=0; j < density_vector.size(); j++)
      if (t_min == 0)
	t_min = density_vector[j];
      else 
	t_min = Minimum(t_min, density_vector[j]);
    
    
    /*
      density_subvector := [density_vector[i] : i in [1..#density_vector] |  density_vector[i] gt 0];
      t_min, t_location := Minimum(density_subvector);
    */
    
    //    /*
    cout << " For the square-class " << t << "Z_" << p << " we have the local density vector ";   
    cout << density_vector << endl;  
    cout << " The minimal value was " << t_min << endl;
      //" at location " << t_location << " of the non-zero subvector ";
    //    */

    t_min.canonicalize();
    const_list[i] = t_min;
    
    
  }
  
  
  /*
    cout <<  "Q = " << Q << "  p = " << p << "  list of squareclasses = " << sqclasslist << endl;
  */
  
    
  //  return const_list, sqclasslist;
  return const_list;
}



/////////////////  The next 3 routines compute numerical things used in Lambda4test


//////////////////////////////////////////////////////////////////////////////
// Secret function to compute the small adjustment factors for Lambda4_test //
//////////////////////////////////////////////////////////////////////////////

mpq_class adjustment_secret(mpz_class m, mpz_class N, mpz_class chi_top) {

  valarray <mpz_class> prime_list;
  prime_list.resize(PrimeDivisors(m).size()); 
  prime_list = PrimeDivisors(m); 

  //  adj_list := [(p-1)/(p+1) : p in PrimeDivisors(m) | (Valuation(N, p) eq 0) and (KroneckerSymbol(chi_top, p) eq -1)];

  mpq_class new_factor;
  new_factor = 1;
  for(size_t i=0; i < prime_list.size(); i++) 
    if ((Valuation(N, prime_list[i]) == 0) && (KroneckerSymbol(chi_top, prime_list[i]) == -1))
      new_factor = new_factor * mpq_class(prime_list[i] - 1, prime_list[i] + 1);
  
  return new_factor;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary function to compute the part of an integer m away from a fixed sequesce S of primes //
///////////////////////////////////////////////////////////////////////////////////////////////////

mpz_class integer_part_away_from_set(mpz_class m, valarray <mpz_class> S) {
  
  mpz_class r;
  r = m;
  
  for(size_t i=0; i<S.size(); i++)
    r = r / (S[i] ^ Valuation(r, S[i]));
  
  return r;
}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute an upper bound for the constant Lambda_4^(hat) numerically using the first few terms of the Eisenstein series //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

mpq_class Lambda4_test(Matrix_mpz Q, valarray <mpq_class> EE) {
  
  mpz_class N;
  N = Q.QFLevel();
  
  cout << " The level is " << N << endl;

  valarray <mpz_class> S;
  S.resize(AnisotropicPrimes(Q).size());
  S = AnisotropicPrimes(Q);
 
  cout << " The anisotropic primes are " << S << endl;
 
  mpz_class chi_top;
  chi_top = CoreDiscriminant(Q.Determinant());
  
  cout << " The character is given by " << chi_top << endl;

  
  valarray <mpq_class> eisratios;
  eisratios.resize(EE.size());
  eisratios[0] = 496;  // This is just a dummy entry...
  for(size_t i=1; i< EE.size(); i++) {
    //    cout << " Starting i = " << i << endl;
    eisratios[i] = EE[i] / (mpq_class(integer_part_away_from_set(i,S)) * adjustment_secret(i, N, chi_top));
  }  

  //    cout << " The Eistein ratios are given by " << eisratios << endl;



  // Return the minimum non-zero entry of Eisratios.
  mpq_class eis_min;
  size_t j_min;

  for(size_t j=1; j < eisratios.size(); j++)
    if (eisratios[j] > 0) {
      eis_min = eisratios[j];
      j_min = j;
    }

  for(size_t j=1; j < eisratios.size(); j++)
    if ((eisratios[j] > 0) && (eisratios[j] < eis_min)) {
      eis_min = eisratios[j];
      j_min = j;
    }

  cout << " The minimum was found at the Eisenstein coefficient m = " << j_min << endl;

  return eis_min;
}





//////////////////////////////////////////////////////////////////// 
// Revised Code for finding bounds for quaternary quadratic forms //
////////////////////////////////////////////////////////////////////

mpq_class new_C4_rational(Matrix_mpz Q) {

  mpz_class N;
  N = Q.QFLevel();
  //print "Level = ", N;

  valarray <mpz_class> N_primes;
  N_primes.resize(PrimeDivisors(N).size());
  N_primes = PrimeDivisors(N);

  //  N_primes_minima := [ -999/2 : i in [1..#N_primes]];


// For each prime and each primitive square class t(Z_p)^2 in Z_p 
// we compute the local factors up to the minimal stable number T', 
// and the constant C_p(T') (together with the adjustment at 
// anisotropic primes).

  mpz_class p;
  valarray <mpq_class> N_primes_minima(N_primes.size());

  for (size_t k=0; k < N_primes.size(); k++){
    p = N_primes[k];
    cout << " p = " << p << endl;

    valarray <mpq_class> temp_C4(PrimSqRepns(p).size());
    temp_C4 = C4_squareclass_constants(Q, p, PrimSqRepns(p));    
    cout << "  and  C4 = " <<  temp_C4;

    N_primes_minima[k] = Minimum(temp_C4);    
    cout << " which gives a minimum of " << N_primes_minima[k] << endl << endl;
  }


  cout << " Primes dividing the level: " <<  N_primes << endl;
  cout << " Associated minima at p: " << N_primes_minima << endl;

  mpz_class d, f;
  mpq_class big_denom, Bigfactor, Finalfactor;

  d = Q.Determinant();
  f = CoreDiscriminant(d);

  big_denom = 1;
  for (size_t k=0; k < N_primes.size(); k++){
    p = N_primes[k];
    big_denom *= (1 - mpq_class(KroneckerSymbol(f,p), p*p));
  }

  Bigfactor = mpq_class(f) * RationalSqrt(mpq_class(f,d)) / (QuadraticBernoulliNumber(2, f.get_si()) * big_denom);


  cout << " Level = " << N << endl;
  cout << " Level primes = " << N_primes << endl;
  cout << " d = " << d << endl;
  cout << " f = " << f << endl;
  cout << " big_denom = " << big_denom << endl;
  cout << " Bigfactor = " << Bigfactor << endl;

  //cout << " Bigfactor = " << Bigfactor << endl;


//print "Bigfactor = ", Bigfactor;
  
  Finalfactor = Bigfactor;
  for (size_t k=0; k < N_primes_minima.size(); k++)
    Finalfactor *= N_primes_minima[k];

  Finalfactor.canonicalize();
  
  return Finalfactor;
}


