////////////////////////////////////////////////////////////////////////////////////////
// This is code to compute the first 10,000 Eisenstein series 
// coefficients of the 6560 quadratic forms in our list...
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



// Make the Eisenstein series

  Q2 := QQ ! Matrix(new_n, new_Q);

  L2 := LatticeWithGram(Q2);

  Gen := GenusRepresentatives(L2);

  AvgTheta := &+[ThetaRing ! ThetaSeries(L, 2*new_precision) / #AutomorphismGroup(L) : L in Gen];

  AvgAuto := &+[ 1 / #AutomorphismGroup(L) : L in Gen];

  Eis2 := AvgTheta / AvgAuto;


// Substitute x for x^2 everywhere

  Eis_half := ThetaRing ! [ Coefficient(Eis2, 2*i) : i in [0..Floor((AbsolutePrecision(Eis2) - 1)/2)]];


// Write the ourput file

  Write(output_filename, Eis_half, "Default");


quit;

