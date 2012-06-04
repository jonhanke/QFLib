
////////////////////////////////////////////////////////////////////////
/// This takes the given (positive definite) quadatic form and reduces
/// it so the diagonal elements are in increasing order, and the
/// off-diagonal entries satisfy 2*|a_ij| <= a_ii for all i < j.
////////////////////////////////////////////////////////////////////////
void Matrix_mpz::ReduceOffDiagonal() {

  // WARNING: This assumes that the input matrix is positive definite! =)

  // Local Variables
  long n = (*this).NumRows();


  // Some basic tests...
  assert((*this).IsSymmetric() == true);
  assert(n >= 1);

  // If we have a 1-variable form, there's nothing to do! =)
  if (n == 1)
    return;


  // Main Loop:
  // ----------
  bool DoneFlag = false;
  while (DoneFlag == false) {

    // First arrange the diagonal elements in increasing order
    while (DoneFlag == false) {
      DoneFlag = true;
      for (long i=1; i<=n-1; i++)
	for (long j=i+1; j<=n; j++)
	  if ((*this)(i,i) > (*this)(j,j)) {
	    (*this).SwapSymmetric(i,j);
	    DoneFlag = false;
	  }
    }
    
    // Run through the row/column of each diagonal entry 
    // to ensure that the basic inequality is satisfied.
    long i=1;
    while ((i <= n-1) && (DoneFlag == true)) {

      long j = i+1;
      while (j <= n) {

	/*
	// DIAGNOSTIC
	cout << " i = " << i << ",  j = " << j << endl;
	cout << (*this) << endl;
	*/


	// Check if the (i,j) entry fails
	if ((*this)(i,i) < abs(2 * (*this)(i,j))) {

	  // Indicate that we must repeat the entire computation again to be sure we're done!
	  DoneFlag = false;

	  // Find the correct multiple of the i-th row/column to add 
	  mpz_class c = - (*this)(i,j) / (*this)(i,i);     // Note: '/' rounds towards zero! =)
	       
	  (*this).AddSymmetric(c, j, i);

	  /*
	  // DIAGNOSTIC
	  cout << " Using c = " << c << endl;
	  cout << (*this) << endl;
	  cout << endl << endl;
	  */
	  
	}

	j++;
      }

      i++;
    }

  }


  // Finally, arrange for the off-diagonal entries to be positive
  for (long i=1; i<=n-1; i++)
    for (long j=i+1; j<=n; j++)
      if ((*this)(i,j) < 0) 
	(*this).MultiplySymmetric(-1,j);
  

  
  return;

}
