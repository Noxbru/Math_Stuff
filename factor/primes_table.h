#ifndef PRIMES_TABLE_H
#define PRIMES_TABLE_H

#include <stdio.h>

unsigned int primes_table_size = 0;
unsigned char *primes_table = NULL;

void generate_primes_table(unsigned int size);

unsigned int get_prime(unsigned int n);

#endif /* end of include guard: PRIMES_TABLE_H */
