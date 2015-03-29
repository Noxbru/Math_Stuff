#ifndef PRIMES_TABLE_H
#define PRIMES_TABLE_H

#include <stdio.h>

extern unsigned int primes_table_size;
#if PRIMES_TABLE_FAT
extern unsigned int *primes_table;
#else
extern unsigned char *primes_table;
#endif

void generate_primes_table(unsigned int size);

unsigned int get_prime(unsigned int n);

#endif /* end of include guard: PRIMES_TABLE_H */
