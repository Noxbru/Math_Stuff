#include <gmp.h>
#include <math.h>

#include <stdlib.h>

#include "factor.h"
#include "prime.h"

/* Please, note that bound_1 and bound_2 are the NUMBER OF PRIMES
 * tried, not that we try primes less than bound_1 for phase 1 and
 * less than bound_2 for phase 2 */
void pollard_p_1(mpz_ptr out, mpz_ptr n, unsigned int bound_1)
{
    const unsigned int bound_2 = 10 * bound_1;
    unsigned int i, size, prime;
    mpz_t aux0, aux1, aux2;
    mpz_t base;

    mpz_inits(aux0, aux1, aux2, NULL);
    mpz_init_set_ui(base, 2);

    generate_primes_table(bound_1);

    mpz_sqrt(aux0, n);
    size = mpz_sizeinbase(aux0, 2);

    i = 1;
    do
    {
        /* top = log_prime(sqrt(n)) */
        prime = get_prime(i);
        unsigned int top = size / log2(prime);
        mpz_ui_pow_ui(aux0, prime, top);
        mpz_powm(aux2, base, aux0, n);

        mpz_swap(base, aux2);

        mpz_sub_ui(aux1, base, 1);

        mpz_gcd(out, aux1, n);

        if(mpz_cmp_ui(out,1) && mpz_cmp(out, n))
            goto out_phase_1;

        i++;
    }
    while (i < bound_1);

    /* Start of Phase 2 of the algorithm */
    mpz_t *aux_base;
    aux_base = malloc(16 * sizeof(mpz_t));

    /* Base will be the 'last base' so here we compute the factors
     * that will take us from one prime to the next */
    mpz_init_set_ui(aux_base[0], 1);
    for(i = 1; i < 16; i++)
    {
        mpz_init(aux_base[i]);
        mpz_powm_ui(aux_base[i], base, 2 * i, n);
    }

    prime = get_prime(bound_1);
    mpz_powm_ui(aux2, base, prime, n);
    mpz_swap(base, aux2);
    mpz_sub_ui(aux1, base, 1);
    mpz_gcd(out, aux1, n);

    /* We might get lucky, although it's very unlikely */
    if(mpz_cmp_ui(out,1) && mpz_cmp(out, n))
        goto out_phase_2;

    /* aux0 = base^p_n */
    mpz_set(aux0, base);
    for(i = prime + 2; bound_1 < bound_2; i+=2)
    {
        if(is_prime_rabin_miller_uint32(i) == PRIME_DEFINITELY_YES)
        {
            unsigned int prime_diff = i - prime;

            while(prime_diff > 30)
            {
                mpz_mul(aux1, aux0, aux_base[15]);
                mpz_mod(aux0, aux1, n);
                prime_diff -= 30;
            }

            mpz_mul(aux1, aux0, aux_base[prime_diff >> 1]);
            mpz_mod(aux0, aux1, n);
            mpz_sub_ui(aux1, aux0, 1);

            mpz_gcd(out, aux1, n);

            if(mpz_cmp_ui(out,1) && mpz_cmp(out, n))
                goto out_phase_2;

            prime = i;
            bound_1++;
        }
    }

out_phase_2:
    for(i = 0; i < 16; i++)
        mpz_clear(aux_base[i]);

out_phase_1:
    mpz_clears(base, aux0, aux1, aux2, NULL);
}
