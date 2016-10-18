#include <gmp.h>

#include "elliptic.h"

int elliptic_curve_sum_montgomery_affine(elliptic_point *p_out,
        elliptic_point *p_in1,
        elliptic_point *p_in2,
        elliptic_context *ctx)
{
    int err = 0;

    if(mpz_cmp(p_in1->x, p_in2->x) == 0)
    {
        if(mpz_cmp(p_in1->y, p_in2->y) == 0)
        {
            return  elliptic_curve_double_montgomery_affine(p_out, p_in1, ctx);
        }
        else
        {
            // set the output point as infinity
            return 0;
        }
    }
    else
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

        mpz_sub(aux1, p_in2->x, p_in1->x);
        mpz_gcdext(lambda, aux2, NULL, aux1, ctx->m);
        if(mpz_cmp_ui(lambda, 1) != 0)
        {
            mpz_set(p_out->x, lambda);
            err = 1;
            goto clean;
        }
        else if(mpz_cmp_ui(aux2, 0) < 0)
            mpz_add(aux2, aux2, ctx->m);

        mpz_sub(aux1, p_in2->y, p_in1->y);
        mpz_mul(lambda, aux1, aux2);        /* lambda = (y2 - y1) / (x2 - x1) */

        mpz_mul(aux1, lambda, lambda);
        mpz_mul(aux2, aux1, ctx->B);
        mpz_sub(aux2, aux2, ctx->A);
        mpz_sub(aux2, aux2, p_in1->x);
        mpz_sub(aux2, aux2, p_in2->x);
        mpz_mod(aux1, aux2, ctx->m);        /* x' = B * lambda^2 - A - x1 - x2 */

        mpz_neg(p_out->y, p_in1->y);
        mpz_sub(aux2, aux1, p_in1->x);
        mpz_submul(p_out->y, aux2, lambda); /* y' = -y1 - lambda * (x' - x1) */

        mpz_swap(p_out->x, aux1);
        mpz_mod(p_out->y, p_out->y, ctx->m);

clean:
        ;
#if !FAT_OBJECTS
        mpz_clears(aux1, aux2, lambda, NULL);
#endif
    }

    return err;
}

int elliptic_curve_double_montgomery_affine(elliptic_point *p_out,
        elliptic_point *p_in,
        elliptic_context *ctx)
{
    int err = 0;

#if FAT_OBJECTS
    mpz_ptr aux1 = ctx->aux1;
    mpz_ptr aux2 = ctx->aux2;
    mpz_ptr lambda = ctx->lambda;
    mpz_ptr nu = ctx->nu;
#else
    mpz_t aux1, aux2;
    mpz_t lambda, nu;
#endif

    if(mpz_cmp_ui(p_in->y, 0) == 0)
    {
        // set the output point as infinity
        return 0;
    }

#if !FAT_OBJECTS
    mpz_inits(aux1, aux2, NULL);
    mpz_inits(lambda, nu, NULL);
#endif

    mpz_add(aux1, p_in->y, p_in->y);
    mpz_mul(aux2, aux1, ctx->B);
    mpz_gcdext(lambda, aux1, NULL, aux2, ctx->m);
    if(mpz_cmp_ui(lambda, 1) != 0)
    {
        mpz_set(p_out->x, lambda);
        err = 1;
        goto clean;
    }
    else if(mpz_cmp_ui(aux1, 0) < 0)
        mpz_add(aux1, aux1, ctx->m);

    mpz_mul(aux2, p_in->x, p_in->x);
    mpz_mul_ui(nu, aux2, 3);
    mpz_mul_2exp(aux2, ctx->A, 1);
    mpz_addmul(nu, aux2, p_in->x);
    mpz_add_ui(nu, nu, 1);
    mpz_mul(lambda, nu, aux1);          /* lambda = (3x^2 +2Ax + 1) / 2By */

    mpz_mul(aux1, lambda, lambda);
    mpz_mul(aux2, aux1, ctx->B);
    mpz_sub(aux2, aux2, ctx->A);
    mpz_submul_ui(aux2, p_in->x, 2);
    mpz_mod(aux1, aux2, ctx->m);        /* x' = B * lambda^2 - A - 2x */

    mpz_neg(p_out->y, p_in->y);
    mpz_sub(aux2, aux1, p_in->x);
    mpz_submul(p_out->y, aux2, lambda); /* y' = -y - lambda *(x' - x) */

    mpz_swap(p_out->x, aux1);
    mpz_mod(p_out->y, p_out->y, ctx->m);

clean:
#if !FAT_OBJECTS
    mpz_clears(aux1, aux2, lambda, nu, NULL);
#endif

    return err;
}

int elliptic_curve_mul_montgomery_affine(elliptic_point *p_out,
        elliptic_point *p_in,
        unsigned int times,
        elliptic_context *ctx)
{
    int err;
    elliptic_point point;
    elliptic_point_init_set(&point, p_in);

    while(times % 2 == 0)
    {
        err = elliptic_curve_double_montgomery_affine(&point, &point, ctx);
        times >>= 1;

        if(err != 0) goto found;
    }

    elliptic_point_set(p_out, &point);

    err = elliptic_curve_double_montgomery_affine(&point, &point, ctx);
    times >>= 1;

    if(err != 0) goto found;

    while(times)
    {
        if(times & 1)
        {
            err = elliptic_curve_sum_montgomery_affine(p_out, p_out,
                    &point, ctx);
            if(err != 0) goto clean;
        }

        err = elliptic_curve_double_montgomery_affine(&point, &point, ctx);
        times >>= 1;
        if(err != 0) goto found;
    }

clean:
    elliptic_point_clear(&point);
    return err;

found:
    mpz_set(p_out->x, point.x);
    goto clean;
}
