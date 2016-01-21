#include <stdio.h>
#include <gmp.h>

#include "elliptic.h"

int main(int argc, const char *argv[])
{
    unsigned int times;
    elliptic_context ctx;
    elliptic_point point;

    if(argc != 7)
        return 1;

    elliptic_curve_init(&ctx);

    mpz_set_str(ctx.A, argv[1], 0);
    mpz_set_str(ctx.B, argv[2], 0);
    mpz_set_str(ctx.m, argv[3], 0);
    ctx.type = ELLIPTIC_CURVE_WEIERSTRASS;

    mpz_init_set_str(point.x, argv[4], 0);
    mpz_init_set_str(point.y, argv[5], 0);
    mpz_init(point.z);
    sscanf(argv[6],"%u",&times);
    point.type = POINT_AFFINE;

    elliptic_curve_mul(&point, &point, times, &ctx);

    gmp_printf("%Zd\t%Zd\n", point.x, point.y);

    elliptic_curve_clear(&ctx);
    elliptic_point_clear(&point);

    return 0;
}
