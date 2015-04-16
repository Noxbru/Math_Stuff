#include <gmp.h>
#include <math.h>

#include "factor.h"
#include "primes_table.h"

void pollard_p_1(mpz_ptr out, mpz_ptr n, unsigned int b)
{
    unsigned int i, j, prime;
    mpz_t aux0, aux1, aux2;
    mpz_t base;

    mpz_inits(aux0, aux1, aux2, NULL);
    mpz_init_set_ui(base, 2);

    generate_primes_table(b);

    mpz_sqrt(aux0, n);

    i = 1;
    do
    {
        /* top = log_prime(sqrt(n)) */
        prime = get_prime(i);
        unsigned int top = mpz_sizeinbase(aux0, 2) /
                           log2(prime);

        mpz_powm_ui(aux1, base, prime, n);

        j = 1;
        do
        {
            mpz_powm_ui(aux2, aux1, prime, n);
            mpz_powm_ui(aux1, aux2, prime, n);
            j+=2;
        }
        while (j < top);

        mpz_swap(base, aux1);

        mpz_sub_ui(aux1, base, 1);

        mpz_gcd(out, aux1, n);

        if(mpz_cmp_ui(out,1) && mpz_cmp(out, n))
            return;

        i++;
    }
    while (i < b);
}
