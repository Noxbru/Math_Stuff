#include <gmp.h>

#include "elliptic.h"

int elliptic_curve_sum_weierstrass_projective(elliptic_point *p_out,
        elliptic_point *p_in1,
        elliptic_point *p_in2,
        elliptic_context *ctx)
{
    if(mpz_cmp(p_in1->x, p_in2->x) == 0 &&
       mpz_cmp(p_in1->y, p_in2->y) == 0 &&
       mpz_cmp(p_in1->z, p_in2->z) == 0)
        return elliptic_curve_double_weierstrass_projective(p_out, p_in1, ctx);

#if FAT_OBJECTS
    mpz_ptr aux1 = ctx->aux1; mpz_ptr aux2 = ctx->aux2;
    mpz_ptr lambda = ctx->lambda;
#else
    mpz_t aux1, aux2;
    mpz_t lambda;

    mpz_inits(aux1, aux2, NULL);
    mpz_inits(lambda, NULL);
#endif
    mpz_t x1z2, y1z2, z1z2, u, uu, v, vv, vvv, a, r;
    mpz_inits(x1z2, y1z2, z1z2, u, uu, v, vv, vvv, a, r, NULL);

    mpz_mul(x1z2, p_in1->x, p_in2->z);
    mpz_mul(y1z2, p_in1->y, p_in2->z);
    mpz_mul(z1z2, p_in1->z, p_in2->z);

    mpz_mul(aux1, p_in1->z, p_in2->y);
    mpz_sub(u, aux1, y1z2);
    mpz_mul(uu, u, u);

    mpz_mul(aux1, p_in1->z, p_in2->x);
    mpz_sub(v, aux1, x1z2);
    mpz_mul(vv, v, v);
    mpz_mul(vvv, vv, v);

    mpz_mul(r, vv, x1z2);

    mpz_mul(aux1, uu, z1z2);
    mpz_sub(a, aux1, vvv);
    mpz_submul_ui(a, r, 2);

    mpz_mul(p_out->x, a, v);
    mpz_mul(p_out->z, vvv, z1z2);

    mpz_sub(aux1, r, a);
    mpz_mul(aux2, aux1, u);
    mpz_mul(lambda, vvv, y1z2);
    mpz_sub(p_out->y, aux2, lambda);

    mpz_mod(p_out->x, p_out->x, ctx->m);
    mpz_mod(p_out->y, p_out->y, ctx->m);
    mpz_mod(p_out->z, p_out->z, ctx->m);

    mpz_clears(x1z2, y1z2, z1z2, u, uu, v, vv, vvv, a, r, NULL);
    return 0;
}

int elliptic_curve_double_weierstrass_projective(elliptic_point *p_out,
        elliptic_point *p_in,
        elliptic_context *ctx)
{
#if FAT_OBJECTS
    mpz_ptr aux1 = ctx->aux1;
    mpz_ptr aux2 = ctx->aux2;
    mpz_ptr lambda = ctx->lambda;
#else
    mpz_t aux1, aux2;
    mpz_t lambda;

    mpz_inits(aux1, aux2, NULL);
    mpz_inits(lambda, NULL);
#endif

    mpz_t xx, zz, w, s, rr, b;
    mpz_inits(xx, zz, w, s, rr, b, NULL);

    mpz_mul(xx, p_in->x, p_in->x);
    mpz_mul(zz, p_in->z, p_in->z);

    mpz_mul(w, ctx->A, zz);
    mpz_addmul_ui(w, xx, 3);

    mpz_mul(aux1, p_in->y, p_in->z);
    mpz_add(s, aux1, aux1);

    mpz_mul(aux2, s, s);
    mpz_mul(p_out->z, s, aux2);

    mpz_mul(aux1, p_in->y, s);
    mpz_mul(rr, aux1, aux1);
    mpz_add(aux1, aux1, p_in->x);
    mpz_mul(b, aux1, aux1);

    mpz_sub(aux1, b, xx);
    mpz_sub(b, aux1, rr);

    mpz_mul(lambda, w, w);
    mpz_submul_ui(lambda, b, 2);

    mpz_mul(p_out->x, lambda, s);

    mpz_sub(aux1, b, lambda);
    mpz_mul(p_out->y, aux1, w);
    mpz_submul_ui(p_out->y, rr, 2);

    mpz_mod(p_out->x, p_out->x, ctx->m);
    mpz_mod(p_out->y, p_out->y, ctx->m);
    mpz_mod(p_out->z, p_out->z, ctx->m);

    mpz_clears(xx, zz, w, s, rr, b, NULL);
    return 0;
}

int elliptic_curve_mul_weierstrass_projective(elliptic_point *p_out,
        elliptic_point *p_in,
        unsigned int times,
        elliptic_context *ctx)
{
    elliptic_point point;
    elliptic_point_init_set(&point, p_in);

    while(times % 2 == 0)
    {
        elliptic_curve_double_weierstrass_projective(&point, &point, ctx);
        times >>= 1;
    }

    elliptic_point_set(p_out, &point);

    elliptic_curve_double_weierstrass_projective(&point, &point, ctx);
    times >>= 1;

    while(times)
    {
        if(times & 1)
            elliptic_curve_sum_weierstrass_projective(p_out, p_out, &point, ctx);

        elliptic_curve_double_weierstrass_projective(&point, &point, ctx);
        times >>= 1;
    }

    return 0;
}
