
// This tests the routines in the old vec_mat.c file:
// --------------------------------------------------

void MatrixTest() {

  /*
  Vector_mpz vv(3);
  cout << " The length of vv is " << vv.length() << endl;
  vv.PrintV(cout);
  */

  //cout << vv << endl;  // Why does this line give so many errors??? =(


  valarray<mpz_class> v;
  v.resize(5);
  v = mpz_class(0);
  cout << " The length of vv is " << v.size() << endl;
  PrintOut(v, v.size());
  

  cout << "\n Now initialize a 3 x 5 matrix: " << endl;
  Matrix_mpz mm(3,5);
  mm.PrintM(cout);

  cout << "\n and take its transpose: " << endl;
  mm.Transpose();
  cout << " Finshed computing the transpose! " << endl;
  mm.PrintM(cout);

  //  cout << vv << endl;


  v[0] = mpz_class(1);
  v[1] = mpz_class(2);
  v[2] = mpz_class(3);
  v[3] = mpz_class(7);

  mm.DiagonalMatrix(v);
  cout << "\n Making a diagonal matrix: " << endl;
  mm.PrintM(cout);
  


  valarray<size_t> I(3);
  I[0] = 1;
  I[1] = 3;
  I[2] = 4;

  //cout << " Wrote the indexing vector" << endl;
  Matrix_mpz mmsub;
  //cout << " Defined the smaller submatrix" << endl;
  mmsub = mm.ExtractSquareSubmatrix(I);
  cout << "\n Extracted the square matrix" << endl;
  mmsub.PrintM(cout);


  size_t i,j;
  cout << "\n Outputting entries with mmsub[i] notation: " << endl;
  for (i=0; i<mmsub.Length(); i++)
    cout << mmsub[i] << " ";
  cout << endl;

  cout << "\n Outputting entries with mmsub(i) notation: " << endl;
  for (i=1; i<=mmsub.Length(); i++)
    cout << mmsub(i) << " ";
  cout << endl;
 

  cout << "\n Outputting entries with mmsub(i,j) notation: " << endl;
  for (i=1; i<=mmsub.NumRows(); i++)
    for (j=1; j<=mmsub.NumCols(); j++)
      cout << mmsub(i,j) << " ";
  cout << endl;



  cout << "\n Changing some entries using mmsub(i,j) notation: " << endl;
  //  mmsub(1,1) = mpz_class(1);  // But we can also use just an integer type! =) 
  mmsub(1,1) = 1;
  mmsub(1,2) = 2;
  mmsub(1,3) = 3;
  mmsub(2,1) = 4;
  mmsub(2,2) = 5;
  mmsub(2,3) = 6;
  mmsub(3,1) = 7;
  mmsub(3,2) = 8;
  mmsub(3,3) = 99;

  mmsub.PrintM(cout);  


  cout << "\n Extracting the [1, 3] submatrix:" << endl; 
  valarray<size_t> J(2);
  J[0] = 1;
  J[1] = 3;
  
  mmsub.ExtractSquareSubmatrix(J).PrintM(cout);
  
  valarray<size_t> Old(3);
  Old[0] = 1;
  Old[1] = 2;
  Old[2] = 3;

  cout << "\n The old index vector ";
  PrintV(Old);
  cout << " is re-indexed to give the new index vector ";
  PrintV(ReindexVectorFromExtraction(Old, J));


  cout << "\n Now let's extract its columns and rows:" << endl;
  for (i=1; i<=mmsub.NumRows(); i++) {
    cout << "  Column " << i << ":  "; 
    PrintV(mmsub.ExtractColumn(i));
  }
  for (j=1; j<=mmsub.NumCols(); j++) {
    cout << "  Row " << j << ":  "; 
    PrintV(mmsub.ExtractRow(j));
  }

}
