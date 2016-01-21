#include <gmp.h>
#include <math.h>

#include "factor.h"
#include "primes_table.h"

#if USE_INTERNAL_MONTGOMERY_MULTIPLICATION
#include "mont.h"
#endif

void pollard_p_1(mpz_ptr out, mpz_ptr n, unsigned int b)
{
    unsigned int i, size, prime;
    mpz_t aux0, aux1, aux2;
    mpz_t base;

#if USE_INTERNAL_MONTGOMERY_MULTIPLICATION
    mont_ctx mctx;
    mont_init(&mctx, n);
#endif

    mpz_inits(aux0, aux1, aux2, NULL);
    mpz_init_set_ui(base, 2);

#if USE_INTERNAL_MONTGOMERY_MULTIPLICATION
    mont_transform(base, base, &mctx);
#endif

    generate_primes_table(b);

    mpz_sqrt(aux0, n);
    size = mpz_sizeinbase(aux0, 2);

    i = 1;
    do
    {
        /* top = log_prime(sqrt(n)) */
        prime = get_prime(i);
        unsigned int top = size / log2(prime);
        mpz_ui_pow_ui(aux0, prime, top);

#if USE_INTERNAL_MONTGOMERY_MULTIPLICATION
        mont_pow_mpz(aux2, base, aux0, &mctx);
        mpz_swap(base, aux2);
        mont_inv_transform(aux2, aux2, &mctx);
        mpz_sub_ui(aux1, aux2, 1);
#else
        mpz_powm(aux2, base, aux0, n);

        mpz_swap(base, aux2);

        mpz_sub_ui(aux1, base, 1);
#endif

        mpz_gcd(out, aux1, n);

        if(mpz_cmp_ui(out,1) && mpz_cmp(out, n))
            return;

        i++;
    }
    while (i < b);

    mpz_clears(base, aux0, aux1, aux2, NULL);
}
