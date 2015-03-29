#include <gmp.h>

#include "primes_table.h"

unsigned int trial_division(mpz_ptr n)
{
    unsigned int i;
    unsigned int prime;

    if(mpz_even_p(n))
        return 2;

    generate_primes_table(16384);

    prime = primes_table[0];

    i = 1;
    do
    {
        prime += primes_table[i];
        if(mpz_divisible_ui_p(n, prime))
            return prime;

        i++;
    }
    while (i < primes_table_size);

    return 1;
}

unsigned int trial_division_ui(mpz_ptr n)
{
    unsigned int prime;

    prime = trial_division(n);
    if(prime != 1)
        return prime;

    prime = get_prime(primes_table_size);

    /* 4294967291 is the biggest prime that fits in a 32
     * bits word */
    {
        if(mpz_divisible_ui_p(n, prime))
            return prime;

        prime+=2;
    }
    while(prime < 4294967292u);

    return 1;
}
