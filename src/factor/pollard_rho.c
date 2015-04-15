#include <gmp.h>

#include "factor.h"

/* Note that the initialization method for seq_1 and seq_2
 * could be different (in many sites it is just used '2')
 * Also, the sequence could advance in different ways, the
 * only restriction is that you can't sum '0' or '-2' */

void pollard_rho(mpz_ptr out, mpz_ptr n)
{
    mpz_t seq_1, seq_2;
    mpz_t aux0, aux1, aux2;

    mpz_init(seq_1);
    mpz_inits(aux0, aux1, aux2, NULL);

    mpz_sqrt(seq_1, n);
    mpz_init_set(seq_2, seq_1);

    do
    {
        /* Advance the first sequence */
        mpz_mul(aux0, seq_1, seq_1);
        mpz_add_ui(aux0, aux0, 1);
        mpz_mod(seq_1, aux0, n);

        /* Advance the second sequence */
        mpz_mul(aux0, seq_2, seq_2);
        mpz_add_ui(aux0, aux0, 1);
        mpz_mul(aux1, aux0, aux0);
        mpz_add_ui(aux1, aux1, 1);
        mpz_mod(seq_2, aux1, n);

        mpz_sub(aux2, seq_1, seq_2);
        mpz_abs(aux2, aux2);
        mpz_gcd(out, aux2, n);
    }
    while (mpz_cmp_ui(out,1) == 0);

    mpz_clears(seq_1, seq_2, aux0, aux1, aux2, NULL);
}
