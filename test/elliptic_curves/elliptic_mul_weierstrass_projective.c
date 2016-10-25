#include <stdio.h>
#include <gmp.h>

#include "elliptic.h"

int main(int argc, const char *argv[])
{
    unsigned int times;
    elliptic_context ctx;
    elliptic_point point;

    if(argc != 8)
        return 1;

    elliptic_curve_init(&ctx);

    mpz_set_str(ctx.A, argv[1], 0);
    mpz_set_str(ctx.B, argv[2], 0);
    mpz_set_str(ctx.m, argv[3], 0);
    ctx.type = ELLIPTIC_CURVE_WEIERSTRASS;

    mpz_init_set_str(point.x, argv[4], 0);
    mpz_init_set_str(point.y, argv[5], 0);
    mpz_init_set_str(point.z, argv[6], 0);
    sscanf(argv[7],"%u",&times);
    point.type = POINT_PROJECTIVE;

    elliptic_curve_mul(&point, &point, times, &ctx);

    mpz_invert(point.z, point.z, ctx.m);
    mpz_mul(point.x, point.x, point.z);
    mpz_mul(point.y, point.y, point.z);
    mpz_mod(point.x, point.x, ctx.m);
    mpz_mod(point.y, point.y, ctx.m);

    gmp_printf("%Zd\t%Zd\t%Zd\n", point.x, point.y, point.z);

    elliptic_curve_clear(&ctx);
    elliptic_point_clear(&point);

    return 0;
}
