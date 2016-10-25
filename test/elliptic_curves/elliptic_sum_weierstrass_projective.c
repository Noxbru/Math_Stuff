#include <stdio.h>
#include <gmp.h>

#include "elliptic.h"

int main(int argc, const char *argv[])
{
    elliptic_context ctx;
    elliptic_point point1;
    elliptic_point point2;

    if(argc != 10)
        return 1;

    elliptic_curve_init(&ctx);

    mpz_set_str(ctx.A, argv[1], 0);
    mpz_set_str(ctx.B, argv[2], 0);
    mpz_set_str(ctx.m, argv[3], 0);
    ctx.type = ELLIPTIC_CURVE_WEIERSTRASS;

    mpz_init_set_str(point1.x, argv[4], 0);
    mpz_init_set_str(point1.y, argv[5], 0);
    mpz_init_set_str(point1.z, argv[6], 0);
    mpz_init_set_str(point2.x, argv[7], 0);
    mpz_init_set_str(point2.y, argv[8], 0);
    mpz_init_set_str(point2.z, argv[9], 0);
    point1.type = POINT_PROJECTIVE;
    point2.type = POINT_PROJECTIVE;

    elliptic_curve_sum(&point1, &point1, &point2, &ctx);

    mpz_invert(point1.z, point1.z, ctx.m);
    mpz_mul(point1.x, point1.x, point1.z);
    mpz_mul(point1.y, point1.y, point1.z);
    mpz_mod(point1.x, point1.x, ctx.m);
    mpz_mod(point1.y, point1.y, ctx.m);

    gmp_printf("%Zd\t%Zd\t%Zd\n", point1.x, point1.y, point1.z);

    elliptic_curve_clear(&ctx);
    elliptic_point_clear(&point1);
    elliptic_point_clear(&point2);

    return 0;
}
