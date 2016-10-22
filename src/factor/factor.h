#ifndef FACTOR_H
#define FACTOR_H

#include <gmp.h>

/* Trial Division based algorithms */

unsigned int trial_division(mpz_ptr n);
unsigned int trial_division_ui(mpz_ptr n);

void elliptic(mpz_ptr out, mpz_ptr n);

void fermat(mpz_ptr out, mpz_ptr n);

void pollard_rho(mpz_ptr out, mpz_ptr n);
void pollard_rho2(mpz_ptr out, mpz_ptr n);
void pollard_rho3(mpz_ptr out, mpz_ptr n);
void pollard_rho4(mpz_ptr out, mpz_ptr n);

void pollard_p_1(mpz_ptr out, mpz_ptr n, unsigned int b);

void quadratic_sieve(mpz_ptr out, mpz_ptr n);

void shanks(mpz_ptr out, mpz_ptr n);

#endif /* end of include guard: FACTOR_H */
