#include <stdlib.h>
#include <math.h>

#include "primes_table.h"

void generate_primes_table(unsigned int size)
{
    unsigned int i;
    unsigned int prime;
    unsigned int last_prime;

    if(size <= primes_table_size)
        return;

    primes_table = realloc(primes_table, size);
    primes_table[0] = 2;
    primes_table[1] = 1;
    primes_table_size = 2;

    last_prime = primes_table[0];
    for(i = 1; i < primes_table_size; i++)
        last_prime += primes_table[i];

    for(prime = last_prime + 2; primes_table_size < size; )
    {
        unsigned int j = 1;
        unsigned int aux = primes_table[0];
        unsigned int sup_limit = sqrt(prime) + 1;

        do
        {
            aux += primes_table[j];
            j++;

            if((prime % aux) == 0)
                goto not_prime;
        }
        while (aux < sup_limit);

        primes_table[primes_table_size] = prime - last_prime;
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

    generate_primes_table(n);

    prime = primes_table[0];
    for(i = 1; i < n; i++)
        prime+=primes_table[i];

    return prime;
}
