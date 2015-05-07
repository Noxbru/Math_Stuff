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

/* One of Brent improvements to Pollard's Rho:
 * change m GDC for m multiplications modulo n plus
 * one GDC */
void pollard_rho2(mpz_ptr out, mpz_ptr n)
{
    unsigned int i, reps;
    mpz_t seq_1, seq_2;
    mpz_t aux0, aux1, aux2, accumulator;

    mpz_init(seq_1);
    mpz_inits(aux0, aux1, aux2, accumulator, NULL);

    mpz_sqrt(seq_1, n);
    mpz_init_set(seq_2, seq_1);

    reps = mpz_sizeinbase(n, 2);

    do
    {
        i = 0;
        mpz_set_ui(accumulator, 1);
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

            mpz_mul(aux0, accumulator, aux2);
            mpz_mod(accumulator, aux0, n);

            i++;
        }
        while (i < reps);

        mpz_gcd(out, accumulator, n);
    }
    while (mpz_cmp_ui(out,1) == 0);

    mpz_clears(seq_1, seq_2, aux0, aux1, aux2, accumulator, NULL);
}

/* More of Brent's improvements to Pollard's Rho:
 * Only calculate F(x) once and only compare some
 * of the values */
void pollard_rho3(mpz_ptr out, mpz_ptr n)
{
    unsigned int i, j, k;
    unsigned int reps;
    mpz_t seq_1, seq_2;
    mpz_t aux0, aux1;
    mpz_t accumulator;

    mpz_init(seq_1);
    mpz_inits(aux0, aux1, accumulator, NULL);

    mpz_sqrt(seq_1, n);
    mpz_init_set(seq_2, seq_1);

    reps = mpz_sizeinbase(n, 2);

    mpz_mul(aux0, seq_2, seq_2);
    mpz_add_ui(aux0, aux0, 1);
    mpz_mod(seq_2, aux0, n);

    j = 1;
    do
    {
        i = j;
        mpz_set(seq_1, seq_2);

        do
        {
            mpz_mul(aux0, seq_2, seq_2);
            mpz_add_ui(aux0, aux0, 1);
            mpz_mod(seq_2, aux0, n);

            j++;
        } while(j < 3*(i + 1)/2);

        do
        {
            k = 0;
            mpz_set_ui(accumulator, 1);
            do
            {
                mpz_mul(aux0, seq_2, seq_2);
                mpz_add_ui(aux0, aux0, 1);
                mpz_mod(seq_2, aux0, n);

                mpz_sub(aux1, seq_1, seq_2);
                mpz_abs(aux1, aux1);

                mpz_mul(aux0, accumulator, aux1);
                mpz_mod(accumulator, aux0, n);

                k++;
                j++;
            }
            while(k < reps && j < 2*i + 1);

            mpz_gcd(out, accumulator, n);

            if(mpz_cmp_ui(out,1) != 0)
                goto out;

        } while(j < 2*i + 1);
    }
    while (1);

out:
    mpz_clears(seq_1, seq_2, aux0, aux1, accumulator, NULL);
}
