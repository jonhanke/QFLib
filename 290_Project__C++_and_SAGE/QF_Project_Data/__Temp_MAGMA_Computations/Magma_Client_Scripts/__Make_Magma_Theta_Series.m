////////////////////////////////////////////////////////////////////////////////////////
// This is code to compute the first precision number of 
// theta series coefficients of the given quadratic form.
// (Actually its the theta series of 2*Q which is guaranteed to be integral and even.)
/////////////////////////////////////////////////////////////////////////////////////////


// Input:
// ------
// new_n := size of n x n matrix
// new_Q := [a_11, a_12, ... , a_nn]  =  matrix as a length n^2 vector
// new_precision := desired precision for the computed Eisenstein series
//
// output_filename := filename where the results should be sent
/////////////////////////////////////////////////////////////////////////////////


// Declare some variables

  Z := Integers();
  Q := Rationals();
  ZZ := RMatrixSpace(Z, new_n, new_n);
  QQ := RMatrixSpace(Q, new_n, new_n);
 
  ThetaRing<q>:=PowerSeriesRing(Q);



// Make the theta series

  Q2 := QQ ! Matrix(new_n, new_Q);

  L2 := LatticeWithGram(Q2);

  Theta := ThetaRing ! ThetaSeries(L2, 2*new_precision);


// Substitute x for x^2 everywhere

  Theta_half := ThetaRing ! [ Coefficient(Theta, 2*i) : i in [0..Floor((AbsolutePrecision(Theta) - 1)/2)]];


// Write the output file

  Write(output_filename, Theta_half, "Default");


quit;

