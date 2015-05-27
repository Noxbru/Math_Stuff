#include <stdio.h>
#include <gmp.h>

#include "elliptic.h"

int main(int argc, const char *argv[])
{
    unsigned int times;
    elliptic_ctx ctx;
    mpz_t x1, y1;

    if(argc != 7)
        return 1;

    mpz_init_set_str(ctx.A, argv[1], 0);
    mpz_init_set_str(ctx.B, argv[2], 0);
    mpz_init_set_str(ctx.m, argv[3], 0);

    mpz_init_set_str(x1, argv[4], 0);
    mpz_init_set_str(y1, argv[5], 0);
    sscanf(argv[6],"%u",&times);

    elliptic_mul(x1, y1, x1, y1, times, ctx);

    gmp_printf("%Zd\t%Zd\n", x1, y1);

    return 0;
}
