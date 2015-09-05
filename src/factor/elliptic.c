#include <gmp.h>
#include <stdio.h>

#include "factor.h"
#include "elliptic.h"
#include "primes_table.h"

void elliptic(mpz_ptr out, mpz_ptr n)
{
    unsigned int i;
    int a;
    elliptic_context e_ctx;
    elliptic_point point;

    /*
     * Numbers obtained by having 3 numbers (n1..n3)
     * and its three squares (s1..s3) and doing:
     *  x1 = (s1 + s2 + s3) / 3
     *  y1 = n1 * n2 * n3
     *  y^2 = (x + s1 - x1) * (x + s2 - x1) * (x + s3 - x1)
     */
    static int coeffs[][4] =
    {
        { 15,  40,   111,  -110},  // (2, 4, 5)
        { 23,  56,   543,  3458},  // (2, 4, 7)
        { 28,  64, -1008, 10368},  // (2, 4, 8)
        { 26,  70,  -507,   506},  // (2, 5, 7)
        { 31,  80,  -927,  5346},  // (2, 5, 8)
        { 39, 112,  -975, -8750},  // (2, 7, 8)
        { 42, 162, -1323,  7722},  // (3, 6, 9)
        { 30, 140,  -291,  1330},  // (4, 5, 7)
        { 35, 160,  -651,  5510},  // (4, 5, 8)
        { 43, 224,  -603, -3402},  // (4, 7, 8)
        { 46, 280,  -387, -1134}   // (5, 7, 8)
    };

    elliptic_point_init(&point);
    elliptic_curve_init(&e_ctx);
    mpz_set(e_ctx.m, n);
    e_ctx.type = ELLIPTIC_CURVE_WEIERSTRASS;
    generate_primes_table(20000);

    for(a = 0; a < 11; a++)
    {
        elliptic_point_set_si_affine(&point, coeffs[a][0], coeffs[a][1]);
        mpz_set_si(e_ctx.A, coeffs[a][2]);
        mpz_set_si(e_ctx.B, coeffs[a][3]);

        for(i = 1; i < 20000; i++)
        {
            if(elliptic_curve_mul(&point, &point, get_prime(i), &e_ctx))
            {
                mpz_set(out, point.x);
                goto clean;
            }
        }
    }

clean:
    elliptic_point_clear(&point);
    elliptic_curve_clear(&e_ctx);
    return;
}
