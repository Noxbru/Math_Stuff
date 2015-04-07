#include <gmp.h>

#include "primes_table.h"

unsigned int trial_division(mpz_ptr n)
{
    unsigned int i;
    unsigned int prime;

    if(mpz_even_p(n))
        return 2;

    /* In both cases we use at least a table with
     * 16 KB of memory */
#if PRIMES_TABLE_FAT
    generate_primes_table(16384);
#else
    generate_primes_table(16384/4);
#endif

    prime = primes_table[0];

    i = 1;
    do
    {
#if PRIMES_TABLE_FAT
        prime = primes_table[i];
#else
        prime += primes_table[i];
#endif
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
