#include <stdio.h>
#include <gmp.h>

#include "elliptic.h"

int main(int argc, const char *argv[])
{
    elliptic_ctx ctx;
    mpz_t x1, y1;
    mpz_t x2, y2;

    if(argc != 8)
        return 1;

    elliptic_init(&ctx);

    mpz_set_str(ctx.A, argv[1], 0);
    mpz_set_str(ctx.B, argv[2], 0);
    mpz_set_str(ctx.m, argv[3], 0);

    mpz_init_set_str(x1, argv[4], 0);
    mpz_init_set_str(y1, argv[5], 0);
    mpz_init_set_str(x2, argv[6], 0);
    mpz_init_set_str(y2, argv[7], 0);

    elliptic_sum(x1, y1, x1, y1, x2, y2, &ctx);

    gmp_printf("%Zd\t%Zd\n", x1, y1);

    return 0;
}
