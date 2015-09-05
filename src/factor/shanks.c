#include <gmp.h>

#include "factor.h"

void shanks(mpz_ptr out, mpz_ptr n)
{
    unsigned int i;
    mpz_t Pi, Pi1;
    mpz_t Qi, Qi1;
    mpz_t b;
    mpz_t aux1, aux2, aux3;

    mpz_inits(Pi, Pi1, Qi, b, aux1, aux2, aux3, NULL);
    mpz_init_set_ui(Qi1, 1);

    mpz_sqrtrem(Pi, Qi, n);
    mpz_set(aux1, Pi);

    i = 1;
    do
    {
        mpz_add(aux2, aux1, Pi);
        mpz_fdiv_q(b, aux2, Qi);

        mpz_swap(Pi1, Pi);
        mpz_mul(aux2, b, Qi);
        mpz_sub(Pi, aux2, Pi1);

        mpz_sub(aux2, Pi1, Pi);
        mpz_addmul(Qi1, aux2, b);
        mpz_swap(Qi, Qi1);

        i++;
    }
    while(!(mpz_perfect_square_p(Qi) && (i % 2 == 0)));

    mpz_sqrt(Qi1, Qi);

    mpz_sub(aux2, aux1, Pi);
    mpz_fdiv_q(b, aux2, Qi1);

    mpz_addmul(Pi, b, Qi1);

    mpz_mul(aux2, Pi, Pi);
    mpz_sub(aux3, n, aux2);
    mpz_fdiv_q(Qi, aux3, Qi1);

    do
    {
        mpz_add(aux2, aux1, Pi);
        mpz_fdiv_q(b, aux2, Qi);

        mpz_swap(Pi1, Pi);
        mpz_mul(aux2, b, Qi);
        mpz_sub(Pi, aux2, Pi1);

        mpz_sub(aux2, Pi1, Pi);
        mpz_addmul(Qi1, aux2, b);
        mpz_swap(Qi, Qi1);
    }
    while(mpz_cmp(Pi, Pi1));

    mpz_gcd(out, Pi, n);

    mpz_clears(Pi, Pi1, Qi, Qi1, b, aux1, aux2, aux3, NULL);
}
