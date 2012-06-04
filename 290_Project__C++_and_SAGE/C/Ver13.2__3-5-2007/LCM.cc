
#include <gmp.h>
#include <gmpxx.h>


////////////////////////////////////
/// Finds the LCM of two integers //
////////////////////////////////////

mpz_class LCM(mpz_class a, mpz_class b) {

  mpz_t l;
  mpz_init(l);
  mpz_lcm(l, a.get_mpz_t(), b.get_mpz_t());

  // Move the LCM
  mpz_class lcm;
  lcm = mpz_class(l);

  // Clear the mpz_t variable
  mpz_clear(l);

  return lcm;
}
