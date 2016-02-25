#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "primes_table.h"
#include "factor.h"

#define tried_numbers 19

static inline void print_bit64(uint64_t u)
{
    char str[65];
    str[64] = '\0';

    for(uint64_t i = 0, id = 1ul << 63;
            i < 64; i++, id >>= 1)
        str[i] = (u & id) ? '1' : '.';
    fputs(str, stdout);
}

static void sort_internal(uint64_t *b,
        mpz_t *x, mpz_t *y, unsigned int size)
{
    unsigned int i, j, m_index;
    uint64_t m;

    for(i = 0; i < size - 1; i++)
    {
        m = b[i];
        m_index = i;
        for(j = i+1; j < size; j++)
        {
            if(b[j] > m)
            {
                m = b[j];
                m_index = j;
            }
        }

        {
            uint64_t lb = b[i]; b[i] = m; b[m_index] = lb;
            mpz_swap(x[i], x[m_index]);
            mpz_swap(y[i], y[m_index]);
        }
    }
}

static void test_numbers(mpz_ptr out,
        mpz_t *x, mpz_t *y,
        unsigned int *indices, unsigned int size,
        mpz_ptr n)
{
    unsigned int i;
    mpz_t x_accumulator, y_accumulator;
    mpz_t aux;

    mpz_init_set_ui(x_accumulator, 1);
    mpz_init_set_ui(y_accumulator, 1);
    mpz_init(aux);

    for(i = 0; i < size; i++)
    {
        mpz_mul(x_accumulator, x_accumulator, x[indices[i]]);
        mpz_mul(y_accumulator, y_accumulator, y[indices[i]]);
    }

    if(mpz_perfect_square_p(y_accumulator))
    {
        printf("ping\n");
        mpz_sqrt(aux, y_accumulator);
        if(mpz_congruent_p(x_accumulator, aux, n))
            printf("Game Over!\n");

        mpz_add(aux, x_accumulator, aux);

        mpz_gcd(out, aux, n);
        gmp_printf("%Zd\n",out);
    }

    mpz_clears(aux, x_accumulator, y_accumulator, NULL);
}

static void solve_by_magic(mpz_ptr out,
        uint64_t *bits,
        mpz_t *x, mpz_t *y, unsigned int size,
        mpz_ptr n)
{
    unsigned int i, j;
    unsigned int *stack_of_indices;
    unsigned int number_of_indices;
    uint64_t b;

    stack_of_indices = malloc(size * sizeof(unsigned int));

    for(i = 0; i < size; i++)
    {
        b = bits[i];
        stack_of_indices[0] = i;
        number_of_indices = 1;

        for(j = stack_of_indices[number_of_indices - 1] + 1;
                j < size || b != 0;
                j++)
        {
            /* We have just removed one index, continue */
            if(bits[j] > b)
                continue;
            /* Trick to check if we can remove a MSB from b:
             * if the same MSB from b is set in bits[j], then
             * b >> 1 will be less than bits[j] */
            if(bits[j] > (b >> 1))
            {
                /* Put index onto stack */
                b ^= bits[j];
                stack_of_indices[number_of_indices] = j;
                number_of_indices++;
            }
            else
            {
                if(number_of_indices == 1)
                    break;

                /* Pop index from stack */
                number_of_indices--;
                j = stack_of_indices[number_of_indices];
                b ^= bits[j];
            }
        }

        if(b == 0)
        {
            printf("FOUND!\n");
            for(j = 0; j < number_of_indices; j++)
            {
                printf("%u\t", stack_of_indices[j]);
            }
            printf("\n");

            test_numbers(out, x, y, stack_of_indices, number_of_indices, n);

            if(mpz_cmp_ui(out,1) && mpz_cmp(out, n))
                goto out;
        }
    }

out:
    free(stack_of_indices);
}

void quadratic_sieve(mpz_ptr out, mpz_t n)
{
    unsigned int i, j, k;

    uint64_t *bits;
    unsigned int *indices;

    mpz_t *relations_x;
    mpz_t *relations_y;

    unsigned int *prime_base;
    unsigned int prime_base_size = 15;

    mpz_t aux0, aux1, aux2, aux3;

    mpz_inits(aux0, aux1, aux2, aux3, NULL);
    relations_x = malloc(tried_numbers * sizeof(mpz_t));
    relations_y = malloc(tried_numbers * sizeof(mpz_t));

    bits = calloc(tried_numbers, sizeof(uint64_t));
    indices = malloc(tried_numbers * sizeof(unsigned int));

    /*
     * We need more primes than twice the size of the base we are
     * going to use because more or less half of them won't be in
     * the prime base. Hopefully, we calculate enough so we don't
     * need to grow the table every time we need a prime
     */
    generate_primes_table(prime_base_size * 2.1);
    prime_base = malloc(prime_base_size * sizeof(unsigned int));
    prime_base[0] = 2;

    /* i iterates over the primes, j puts them in the base */
    for(i = 0, j = 1; j < prime_base_size; i++)
    {
        unsigned int prime = get_prime(i + 1);
        int kronecker_symbol = mpz_kronecker_ui(n, prime);

        /* huh? there was a little prime as factor */
        if(kronecker_symbol == 0)
        {
            mpz_set_ui(out, prime);
            goto out2;
        }

        /* Prime for the base :) */
        if(kronecker_symbol == 1)
        {
            prime_base[j] = prime;
            j++;
        }
    }

    printf("The prime base is:\t");
    for(i = 0; i < prime_base_size; i++)
    {
        printf("%u\t",prime_base[i]);
    }
    printf("\n");

    unsigned int poly_factor = 1;
    unsigned int numbers_per_poly = 100000;
    for(i = 0; i < tried_numbers; )
    {
        mpz_mul_ui(aux1, n, poly_factor);

        mpz_sqrt(aux0, aux1);
        mpz_add_ui(aux0, aux0, 1);
        mpz_mul(aux1, aux0, aux0);

        mpz_sub(aux2, aux1, n);

        /*gmp_printf("%Zd\n", aux0);*/

        for(j = 0; j < numbers_per_poly && i < tried_numbers; j++)
        {
            /* We need to keep 'aux2' around */
            mpz_set(aux3, aux2);

            for(k = 0; k < prime_base_size; k++)
            {
                while(mpz_divisible_ui_p(aux3, prime_base[k]))
                {
                    mpz_divexact_ui(aux3, aux3, prime_base[k]);
                    bits[i] ^= 1ul << k;
                }
            }

            if(mpz_cmp_ui(aux3, 1) != 0)
                bits[i] = 0;
            /* One usual relation, add it to the list */
            else if(bits[i] != 0)
            {
                mpz_init_set(relations_x[i], aux0);
                mpz_init_set(relations_y[i], aux2);
                i++;
            }
            /* Quadratic relation! nice! */
            else
            {
                mpz_sqrt(aux3, aux2);
                mpz_add(aux3, aux3, aux0);
                mpz_gcd(out, aux3, n);

                if(mpz_cmp_ui(out,1) && mpz_cmp(out, n))
                    goto out;
            }

            mpz_addmul_ui(aux2, aux0, 2);
            mpz_add_ui(aux2, aux2, 1);
            mpz_add_ui(aux0, aux0, 1);
        }

        poly_factor++;
        printf("POLYNOMIAL CHANGE!\n");
    }

    sort_internal(bits, relations_x, relations_y, tried_numbers);

    for(i = 0; i < tried_numbers; i++)
    {
        gmp_printf("%2d:\t%10Zd\t%15Zd\t", i, relations_x[i], relations_y[i]);
        print_bit64(bits[i]);
        printf("\n");
    }

    solve_by_magic(out, bits, relations_x, relations_y, tried_numbers, n);

out:
    for(i = 0; i < tried_numbers; i++)
    {
        mpz_clear(relations_x[i]);
        mpz_clear(relations_y[i]);
    }

out2:
    for(i = 1; i < tried_numbers; i++)
    {
        bits[0] |= bits[i];
    }
    print_bit64(bits[0]); printf("\n");

    free(relations_x);
    free(relations_y);
    free(bits);
    free(indices);
    free(prime_base);
    mpz_clears(aux0, aux1, aux2, aux3, NULL);
}
