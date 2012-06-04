
// This tests all of the routines in vector.cc

void VectorTest() {

  cout << "\n\n Testing the Vector Routines:" << endl;
  cout << " ----------------------------" << endl;

  valarray<size_t> u(4), v(3), w(5);

  u[0] = 2;
  u[1] = 4;
  u[2] = 6;
  u[3] = 7;

  v[0] = 1;
  v[1] = 3;
  v[2] = 5;

  w[0] = 8;
  w[1] = 4;
  w[2] = 2;
  w[3] = 3;
  w[4] = 7;

  // Here are the starting vectors
  cout << "\n Starting with the vectors: " << endl;
  cout << "  u = ";
  PrintV(u);

  cout << "  v = ";
  PrintV(v);

  cout << "  w = ";
  PrintV(w);


  // Tests OrderVector routine
  cout << "\n Order each of the vectors: " << endl;
  cout << "  OrderVector(u):  ";
  PrintV(OrderVector(u));

  cout << "  OrderVector(v):  ";
  PrintV(OrderVector(v));

  cout << "  OrderVector(w):  ";
  PrintV(OrderVector(w));


  // Tests IsDisjointOrdered routine (need the vectors to be ordered!)
  cout << "\n Check which of the 3 vectors are disjoint: " << endl;
  cout << "  u and v?  " << IsDisjointOrdered(u,v) << endl;
  cout << "  u and w?  " << IsDisjointOrdered(u,w) << endl;
  cout << "  v and w?  " << IsDisjointOrdered(v,w) << endl;
  cout << "  v and OrderVector(w)?  " << IsDisjointOrdered(v,OrderVector(w)) << endl;


  // Tests the MakeVector routine
  cout << "\n Made the vector: ";
  PrintV(MakeVector(7,1,2));


  // Tests the VectorAppend Routine
  cout << "\n Appending some vectors:  " << endl;
  cout << "  Append u and v:  ";
  PrintV(VectorAppend(u,v));

  cout << "  Append u and w:  ";
  PrintV(VectorAppend(u,w));

  cout << "  Append v and w:  ";
  PrintV(VectorAppend(v,w));


  // Tests the VectorUnion routine
  cout << "\n Unions of some vectors:  " << endl;
  cout << "  Union of u and v:  ";
  PrintV(VectorUnion(u,v));

  cout << "  Union of u and w:  ";
  PrintV(VectorUnion(u,w));

  cout << "  Union of v and w:  ";
  PrintV(VectorUnion(v,w));


  /*
  // Tests the VectorUnionOrdered routine -- (currently the same as VectorUnion)
  cout << "\n Ordered Unions of some vectors:  " << endl;
  cout << "  Ordered Union of u and v:  ";
  PrintV(VectorUnionOrdered(u,v));

  cout << "  Ordered Union of u and w:  ";
  PrintV(VectorUnionOrdered(u,w));

  cout << "  Ordered Union of v and w:  ";
  PrintV(VectorUnionOrdered(v,w));
  */


  // Tests the VectorIntersection routine
  cout << "\n Intersections of some vectors:  " << endl;
  cout << "  Intersection of u and v:  ";
  PrintV(VectorIntersection(u,v));

  cout << "  Intersection of u and w:  ";
  PrintV(VectorIntersection(u,w));

  cout << "  Intersection of v and w:  ";
  PrintV(VectorIntersection(v,w));


  /*
  // Tests the VectorIntersectionOrdered routine -- (currently the same as VectorIntersection)
  cout << "\n Ordered Intersections of some vectors:  " << endl;
  cout << "  Ordered Intersection of u and v:  ";
  PrintV(VectorIntersectionOrdered(u,v));

  cout << "  Ordered Intersection of u and w:  ";
  PrintV(VectorIntersectionOrdered(u,w));

  cout << "  Ordered Intersection of v and w:  ";
  PrintV(VectorIntersectionOrdered(v,w));
  */


  // Tests the VectorRemove routine
  cout << "\n Remove some vectors from each other:  " << endl;
  cout << "  Remove u from v:  ";
  PrintV(VectorRemove(v,u));

  cout << "  Remove u from w:  ";
  PrintV(VectorRemove(w,u));

  cout << "  Remove v from w:  ";
  PrintV(VectorRemove(w,v));


  // Tests the VectorComplement routine
  cout << "\n Complements of some vectors:  " << endl;
  cout << "  Complement of u by v:  ";
  PrintV(VectorComplement(u,v));

  cout << "  Complement of u by w:  ";
  PrintV(VectorComplement(u,w));

  cout << "  Complement of v by w:  ";
  PrintV(VectorComplement(v,w));

  cout << "  Complement of w by u:  ";
  PrintV(VectorComplement(w,u));

  cout << "  Complement of w by v:  ";
  PrintV(VectorComplement(w,v));


  /*
  // Tests the VectorComplementOrdered routine -- (currently the same as VectorComplement)
  cout << "\n Ordered Complements of some vectors:  " << endl;
  cout << "  Ordered Complement of u by v:  ";
  PrintV(VectorComplementOrdered(u,v));

  cout << "  Ordered Complement of u by w:  ";
  PrintV(VectorComplementOrdered(u,w));

  cout << "  Ordered Complement of v by w:  ";
  PrintV(VectorComplementOrdered(v,w));

  cout << "  Ordered Complement of w by u:  ";
  PrintV(VectorComplementOrdered(w,u));

  cout << "  Ordered Complement of w by v:  ";
  PrintV(VectorComplementOrdered(w,v));
  */


  // Tests the CheckVectorIntersection routine 
  cout << "\n Check which of the 3 vectors intersect: " << endl;
  cout << "  u and v?  " << CheckVectorIntersection(u,v) << endl;
  cout << "  u and w?  " << CheckVectorIntersection(u,w) << endl;
  cout << "  v and w?  " << CheckVectorIntersection(v,w) << endl;


  /*
  // Tests the CheckVectorIntersectionOrdered routine -- (currently the same as CheckVectorIntersection)
  cout << "\n Check which of the 3 ordered vectors intersect: " << endl;
  cout << "  u and v?  " << CheckVectorIntersectionOrdered(u,v) << endl;
  cout << "  u and w?  " << CheckVectorIntersectionOrdered(u,w) << endl;
  cout << "  v and w?  " << CheckVectorIntersectionOrdered(v,w) << endl;
  cout << "  v and OrderVector(w)?  " << CheckVectorIntersectionOrdered(v,OrderVector(w)) << endl;
  */


  // Tests the two GCD routines (pair of numbers, and for a vector of numbers)
  cout << "\n Check the GCDs of a pair of numbers : " << endl;
  cout << "  2 and 7?  " << GCD(2,7) << endl;
  cout << "  -2 and 7?  " << GCD(-2,7) << endl;
  cout << "  -2 and -7?  " << GCD(-2,-7) << endl;
  cout << "  -12 and 9?  " << GCD(-12,9) << endl;
  cout << "  -0 and 3?  " << GCD(-0,3) << endl;
  cout << "  -0 and 0?  " << GCD(-0,0) << endl;

  valarray<mpz_class> u1(4), v1(3), w1(5);

  u1[0] = 24;
  u1[1] = 4;
  u1[2] = -6;
  u1[3] = 18;

  v1[0] = 7;
  v1[1] = 3;
  v1[2] = 5;

  w1[0] = 54;
  w1[1] = 0;
  w1[2] = 81;
  w1[3] = 3;
  w1[4] = -12;

  cout << "\n Check the GCD of a vector of numbers : " << endl;
  cout << "  GCD of ";
  PrintV(u1);
  cout << " is " << GCD(u1) << endl;

  cout << "  GCD of ";
  PrintV(v1);
  cout << " is " << GCD(v1) << endl;

  cout << "  GCD of ";
  PrintV(w1);
  cout << " is " << GCD(w1) << endl;

}

