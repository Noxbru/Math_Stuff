#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include <gmp.h>

typedef struct
{
    mpz_t A, B;
    mpz_t m;
} elliptic_ctx;

int elliptic_sum(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in1, mpz_ptr y_in1,
        mpz_ptr x_in2, mpz_ptr y_in2,
        elliptic_ctx ctx);

int elliptic_double(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in, mpz_ptr y_in,
        elliptic_ctx ctx);

int elliptic_mul(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in, mpz_ptr y_in,
        unsigned int times,
        elliptic_ctx ctx);

#endif /* end of include guard: ELLIPTIC_H */
