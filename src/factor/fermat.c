#include <gmp.h>
#include <stdio.h>

#include "factor.h"

void fermat(mpz_ptr out, mpz_ptr n)
{
    mpz_t x, y, x2, y2;

    mpz_inits(x, y, x2, y2, NULL);

    mpz_sqrt(x, n);

    do
    {
        mpz_add_ui(x, x, 1);
        mpz_mul(x2, x, x);

        mpz_sub(y2, x2, n);

        if(mpz_perfect_square_p(y2))
        {
            mpz_sqrt(y, y2);
            mpz_sub(out, x, y);

            mpz_clears(x, y, x2, y2, NULL);
            return;
        }
    }
    while (1);
}
