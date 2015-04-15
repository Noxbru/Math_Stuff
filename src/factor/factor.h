#ifndef FACTOR_H
#define FACTOR_H

#include <gmp.h>

/* Trial Division based algorithms */

unsigned int trial_division(mpz_ptr n);
unsigned int trial_division_ui(mpz_ptr n);

void fermat(mpz_ptr out, mpz_ptr n);

void pollard_rho(mpz_ptr out, mpz_ptr n);

#endif /* end of include guard: FACTOR_H */
