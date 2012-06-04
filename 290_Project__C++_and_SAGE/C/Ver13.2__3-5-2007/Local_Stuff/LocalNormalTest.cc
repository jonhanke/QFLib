
void LocalNormalTest() {

  // Checks our custom exponentiation operator which
  // allows us to evaluate: mpz_class^(unsigned long)
  mpz_class bb = 3, ans;
  unsigned long pp = 2;
  ans = bb ^ pp;

  cout << "\n Computing " << bb << "^" << pp << endl;
  cout << "  which gives " << (bb^pp) << endl;        // Note: These parentheses are necessary!
  cout << "  which gives " << ans << endl;    


  // Checks the ZeroMatrix routine
  cout << "\n The 2 x 6 zero matrix is: " << endl << ZeroMatrix(2,6) << endl;

  // Checks the IdentityMatrix routine
  cout << "\n The 4 x 4 identity matrix is: " << endl << IdentityMatrix(4) << endl;


  // Make a 3 x 4 matrix M for future testing...
  Matrix_mpz M(3,4);

  M(1,1) = 1;
  M(1,2) = 2;
  M(1,3) = 3;
  M(1,4) = 4;
  M(2,1) = 5;
  M(2,2) = 6;
  M(2,3) = 7;
  M(2,4) = 8;
  M(3,1) = 9;
  M(3,2) = 10;
  M(3,3) = 11;
  M(3,4) = 12;
  cout << "\n Now start with the 3 x 5 matrix M: " << endl << M;


  // Checks the ResizeMatrix routine
  cout << "\n  We resize it to the 2 x 3 matrix: " << endl << ResizeMatrix(M,2,3); // << endl;
  cout << "  and we resize it to the 4 x 6 matrix: " << endl << ResizeMatrix(M,4,6); // << endl;


  // Make a 3 x 3 symmetric matrix S for future testing...
	Matrix_mpz S(3,3);
	// CAN'T GET EITHER OF THESE TO WORK... =(
	//	S = ResizeMatrix(M,3,3) + ResizeMatrix(M,3,3);
	//	S = ResizeMatrix(M,3,3) + ResizeMatrix(M,3,3).Transpose();
  S(1,1) = 1;
  S(1,2) = S(2,1) = 2;
  S(1,3) = S(3,1) = 3;
  S(2,2) = 4;
  S(2,3) = S(3,2) = 5;
  S(3,3) = 6;
  cout << "\n\n  Start with the symmetrix 3 x 3 matrix S: " << endl << S;


  // Checks the SwapSymmetric routine
  cout << "\n  Now Swap the rows/columns 1 and 3 of S: " << endl << SwapSymmetric(S,1,3);

  // Checks the MultiplySymmetric routine
  cout << "\n  Now Multiply the row/column 2 by 4 of S: " << endl << MultiplySymmetric(S,4,2);

  // Checks the AddSymmetric routine
  cout << "\n  Now Add 2 x row/column 1 row/column 3 of S: " << endl << AddSymmetric(S,2,1,3);


  
  // Checks the LocalNormal routine
    cout << "\n  The Local Normal form of S at p=2 is: " << endl << LocalNormal(S,2);
  cout << "\n  The Local Normal form of S at p=3 is: " << endl << LocalNormal(S,3);
  cout << "\n  The Local Normal form of S at p=5 is: " << endl << LocalNormal(S,5);
  cout << "\n  The Local Normal form of S at p=7 is: " << endl << LocalNormal(S,7);
  

  
  // Checks the Determinant member function
  cout << "\n  The matrix S is: \n" << S << endl;
  cout << "\n  The Determinant of S is: " << S.Determinant() << endl;

  Matrix_mpz S1 = ResizeMatrix(S,2,2);
  cout << "\n  The Determinant of \n" << S1 <<"  is: " << S1.Determinant() << endl;

  S1(2,2) = 5;
  cout << "\n  The Determinant of \n" << S1 <<"  is: " << S1.Determinant() << endl;

}


