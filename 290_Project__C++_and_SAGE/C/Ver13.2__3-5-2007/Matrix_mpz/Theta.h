
#if !defined(THETA_H)
#define THETA_H



/*
#include "Matrix_mpz.h"
#include "../New_Boolean_Theta_Class/boolean_theta_class.h"  // This is only needed for our theta function computation...!
*/


  // ==================================== Theta.cc =================================================
  PowerSeries<mpz_class> Theta_PARI_1(const Matrix_mpz & Q, const mpz_class & precision);
  boolean_theta Theta_PARI_2(const Matrix_mpz & Q, const mpz_class & precision);
  PowerSeries<mpz_class> Theta_PARI_1_new(const Matrix_mpz & Q, const mpz_class & precision);
  boolean_theta Theta_PARI_2_new(const Matrix_mpz & Q, const mpz_class & precision);

  boolean_theta Theta_PARI_2_new_Approximate(const Matrix_mpz & Q, const mpz_class & precision, const float & approx_size = 500.);
  boolean_theta Theta_PARI_3_new_Approximate_Ternary(const Matrix_mpz & Q, const mpz_class & precision, const float & approx_size = 500.);
  boolean_theta Theta_PARI_3_new_Approximate_Ternary_mpz(const Matrix_mpz & Q, const mpz_class & precision, const float & approx_size = 500.);

  PowerSeries<mpz_class> Theta1(const Matrix_mpz & Q, const mpz_class & precision);
  boolean_theta Theta2(const Matrix_mpz & Q, const mpz_class & precision);
  mpz_class _Jagy_IntSqrt(const mpz_class & n);





// Theta function computation
//#include "Theta.cc"

#endif 
