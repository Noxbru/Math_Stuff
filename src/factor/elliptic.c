#include <gmp.h>
#include <stdio.h>

#include "factor.h"
#include "elliptic.h"
#include "primes_table.h"

void elliptic(mpz_ptr out, mpz_ptr n)
{
    unsigned int i;
    int a;
    elliptic_ctx e_ctx;
    mpz_t x, y;

    mpz_inits(x, y, e_ctx.A, e_ctx.B, NULL);
    mpz_init_set(e_ctx.m, n);
    generate_primes_table(20000);

    for(a = 2; a < 8; a++)
    {
        mpz_set_ui(x, 1);
        mpz_set_ui(y, 1);
        mpz_set_ui(e_ctx.A, a);
        mpz_set_si(e_ctx.B, -a);

        for(i = 1; i < 20000; i++)
        {
            if(elliptic_mul(x, y, x, y, i, &e_ctx))
            {
                mpz_set(out, x);
                goto clean;
            }
        }
    }

clean:
    mpz_clears(x, y, e_ctx.A, e_ctx.B, e_ctx.m, NULL);
    return;
}
