#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "primes_table.h"

unsigned int primes_table_size = 0;
#if PRIMES_TABLE_FAT
unsigned int *primes_table = NULL;
#else
unsigned char *primes_table = NULL;
#endif

void generate_primes_table(unsigned int size)
{
    unsigned int i;
    unsigned int prime;
    unsigned int last_prime;

    if(size <= primes_table_size)
        return;

#if PRIMES_TABLE_FAT
    primes_table = realloc(primes_table, size * sizeof(unsigned int));
    primes_table[1] = 3;
#else
    primes_table = realloc(primes_table, size * sizeof(unsigned char));
    primes_table[1] = 1;
#endif
    primes_table[0] = 2;
    primes_table_size = 2;

#if PRIMES_TABLE_FAT
    last_prime = primes_table[primes_table_size - 1];
#else
    last_prime = primes_table[0];
    for(i = 1; i < primes_table_size; i++)
        last_prime += primes_table[i];
#endif

    for(prime = last_prime + 2; primes_table_size < size; )
    {
        unsigned int j = 1;
        unsigned int aux = primes_table[0];
        unsigned int sup_limit = sqrt(prime) + 1;

        do
        {
#if PRIMES_TABLE_FAT
            aux = primes_table[j];
#else
            aux += primes_table[j];
#endif
            j++;

            if((prime % aux) == 0)
                goto not_prime;
        }
        while (aux < sup_limit);

#if PRIMES_TABLE_FAT
        primes_table[primes_table_size] = prime;
#else
        primes_table[primes_table_size] = prime - last_prime;
#endif
        primes_table_size++;

        last_prime = prime;

not_prime:
        prime+=2;
    }
}

unsigned int get_prime(unsigned int n)
{
    unsigned int i;
    unsigned int prime;

    assert(n);

    generate_primes_table(n);

#if PRIMES_TABLE_FAT
    return primes_table[n - 1];
#else
    prime = primes_table[0];
    for(i = 1; i < n; i++)
        prime+=primes_table[i];

    return prime;
#endif
}
