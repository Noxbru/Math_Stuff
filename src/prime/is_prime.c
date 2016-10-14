#include <math.h>
#include <inttypes.h>

#include "prime.h"

unsigned int is_prime_table(unsigned int n)
{
    unsigned int last_prime;
    unsigned int i;

    if(n == 2) return PRIME_DEFINITELY_YES;

#if PRIMES_TABLE_FAT
    /* Quick exit if n is greater than last_prime^2 */
    last_prime = primes_table[primes_table_size - 1];

    if(last_prime * last_prime < n)
        return PRIME_UNKNOWN;
#endif

    /* n between 2 and last_prime^2
     * We could do better in the case of PRIMES_TABLE_FAT
     * for 0 < n < last_prime by using bisection,
     * but this function should not be called much, so we
     * prefer to keep more common code*/
    last_prime = 2;
    for(i = 1; i < primes_table_size; i++)
    {
#if PRIMES_TABLE_FAT
        unsigned int next_prime = primes_table[i];
#else
        unsigned int next_prime = last_prime + primes_table[i];
#endif

        if(n == next_prime)      return PRIME_DEFINITELY_YES;
        if(n < next_prime)       return PRIME_DEFINITELY_NO;
        if(n % next_prime == 0)  return PRIME_DEFINITELY_NO;

#if !PRIMES_TABLE_FAT
        last_prime = next_prime;
#endif
    }

#if !PRIMES_TABLE_FAT
    /* Case with n greater than last_prime^2 */
    if(last_prime * last_prime < n)
        return PRIME_UNKNOWN;
#endif

    /* Here we reach in the case last_prime < n < last_prime^2
     * in the case we haven't found a factor among the primes
     * of the table, so the number must be prime */
    return PRIME_DEFINITELY_NO;
}

static uint32_t powm_uint32(uint32_t base, uint32_t exp, uint32_t mod)
{
    uint64_t b = base;
    uint64_t r = 1;

    while(exp)
    {
        if((exp & 1) == 1)
            r *= b;
        b *= b;

        b %= mod;
        r %= mod;
        exp >>= 1;
    }

    return r;
}

unsigned int is_strong_pseudoprime(unsigned int base, unsigned int n)
{
    unsigned int b, t, q;
    unsigned int i;

    q = n - 1;

    t = __builtin_ctz(q);
    q >>=t;

    b = powm_uint32(base, q, n);

    if(b == 1)
        return 1;

    for(i = 1; i < t; i++)
    {
        if(b == 1 || b == n - 1)
            break;

        b = ((uint64_t) b*b) % n;
    }

    return b == n - 1;
}

unsigned int is_prime_rabin_miller_uint32(uint32_t n)
{
    if(n ==   25326001 ||
       n ==  161304001 ||
       n ==  960946321 ||
       n == 1157839381 ||
       n == 3215031751 ||
       n == 3697278427)
        return PRIME_DEFINITELY_NO;

    if(is_strong_pseudoprime(2, n) == 0 ||
       is_strong_pseudoprime(3, n) == 0 ||
       is_strong_pseudoprime(5, n) == 0)
        return PRIME_DEFINITELY_NO;


    return PRIME_DEFINITELY_YES;
}
