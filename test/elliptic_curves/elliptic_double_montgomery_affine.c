#include <stdio.h>
#include <gmp.h>

#include "elliptic.h"

int main(int argc, const char *argv[])
{
    elliptic_context ctx;
    elliptic_point point;

    if(argc != 6)
        return 1;

    elliptic_curve_init(&ctx);

    mpz_set_str(ctx.A, argv[1], 0);
    mpz_set_str(ctx.B, argv[2], 0);
    mpz_set_str(ctx.m, argv[3], 0);
    ctx.type = ELLIPTIC_CURVE_MONTGOMERY;

    mpz_init_set_str(point.x, argv[4], 0);
    mpz_init_set_str(point.y, argv[5], 0);
    mpz_init(point.z);
    point.type = POINT_AFFINE;

    elliptic_curve_double(&point, &point, &ctx);

    gmp_printf("%Zd\t%Zd\n", point.x, point.y);

    elliptic_curve_clear(&ctx);
    elliptic_point_clear(&point);

    return 0;
}
