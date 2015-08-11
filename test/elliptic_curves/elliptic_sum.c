#include <stdio.h>
#include <gmp.h>

#include "elliptic.h"

int main(int argc, const char *argv[])
{
    elliptic_context ctx;
    elliptic_point point1;
    elliptic_point point2;

    if(argc != 8)
        return 1;

    elliptic_init(&ctx);

    mpz_set_str(ctx.A, argv[1], 0);
    mpz_set_str(ctx.B, argv[2], 0);
    mpz_set_str(ctx.m, argv[3], 0);

    mpz_init_set_str(point1.x, argv[4], 0);
    mpz_init_set_str(point1.y, argv[5], 0);
    mpz_init_set_str(point2.x, argv[6], 0);
    mpz_init_set_str(point2.y, argv[7], 0);
    mpz_init(point1.z);
    mpz_init(point2.z);
    point1.type = POINT_AFFINE;
    point2.type = POINT_AFFINE;

    elliptic_sum(&point1, &point1, &point2, &ctx);

    gmp_printf("%Zd\t%Zd\n", point1.x, point1.y);

    return 0;
}
