#include <gmp.h>

#include "factor.h"

void fermat(mpz_ptr out, mpz_ptr n)
{
    mpz_t aux0, aux1, aux2;

    mpz_inits(aux0, aux1, aux2, NULL);

    mpz_sqrt(aux0, n);

    do
    {
        mpz_add_ui(aux0, aux0, 1);
        mpz_mul(aux1, aux0, aux0);

        mpz_sub(aux1, aux1, n);

        if(mpz_perfect_square_p(aux1))
        {
            mpz_sqrt(aux2, aux1);
            mpz_sub(out, aux0, aux2);

            mpz_clears(aux0, aux1, aux2, NULL);
            return;
        }
    }
    while (1);
}
