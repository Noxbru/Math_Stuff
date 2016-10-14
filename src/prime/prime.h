#ifndef PRIMES_TABLE_H
#define PRIMES_TABLE_H

#include <stdio.h>
#include <inttypes.h>

extern unsigned int primes_table_size;
#if PRIMES_TABLE_FAT
extern unsigned int *primes_table;
#else
extern unsigned char *primes_table;
#endif

enum prime_status
{
    PRIME_DEFINITELY_NO,
    PRIME_DEFINITELY_YES,
    PRIME_MAYBE_YES,
    PRIME_UNKNOWN
};

void generate_primes_table(unsigned int size);

unsigned int get_prime(unsigned int n);

unsigned int is_prime_table(unsigned int n);
unsigned int is_prime_rabin_miller_uint32(uint32_t n);

unsigned int is_strong_pseudoprime(unsigned int base, unsigned int n);

#endif /* end of include guard: PRIMES_TABLE_H */
