#include <gmp.h>

#include "elliptic.h"

int elliptic_sum(elliptic_point *p_out,
        elliptic_point *p_in1,
        elliptic_point *p_in2,
        elliptic_ctx *ctx)
{
    int err;

    if(mpz_cmp(p_in1->x, p_in2->x) == 0)
    {
        if(mpz_cmp(p_in1->y, p_in2->y) == 0)
        {
            return elliptic_double(p_out, p_in1, ctx);
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
        mpz_ptr nu = ctx->nu;
#else
        mpz_t aux1, aux2;
        mpz_t lambda, nu;

        mpz_inits(aux1, aux2, NULL);
        mpz_inits(lambda, nu, NULL);
#endif

        mpz_sub(aux1, p_in2->x, p_in1->x);
        err = mpz_invert(aux2, aux1, ctx->m);
        if(err == 0)
        {
            mpz_gcd(p_out->x, aux1, ctx->m);
            goto clean;
        }

        mpz_sub(aux1, p_in2->y, p_in1->y);
        mpz_mul(lambda, aux1, aux2);    /* lambda = (y2 - y1) / (x2 - x1) */

        mpz_mul(aux1, p_in1->x, lambda);
        mpz_sub(nu, p_in1->y, aux1);       /* nu = y1 - lambda * x1 */

        mpz_mul(aux1, lambda, lambda);
        mpz_sub(aux1, aux1, p_in1->x);
        mpz_sub(aux1, aux1, p_in2->x);
        mpz_swap(p_out->x, aux1);          /* x' = lambda^2 - x1 - x2 */

        mpz_neg(p_out->y, nu);
        mpz_submul(p_out->y, p_out->x, lambda);   /* y' = -(lambda * x' + nu) */

        mpz_mod(p_out->x, p_out->x, ctx->m);
        mpz_mod(p_out->y, p_out->y, ctx->m);

clean:
        ;
#if !FAT_OBJECTS
        mpz_clears(aux1, aux2, lambda, nu, NULL);
#endif
    }

    return !err;
}

int elliptic_double(elliptic_point *p_out,
        elliptic_point *p_in,
        elliptic_ctx *ctx)
{
    int err;

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
    err = mpz_invert(aux2, aux1, ctx->m);
    if(err == 0)
    {
        mpz_gcd(p_out->x, aux1, ctx->m);
        goto clean;
    }

    mpz_mul(aux1, p_in->x, p_in->x);
    mpz_mul_ui(aux1, aux1, 3);
    mpz_add(aux1, aux1, ctx->A);
    mpz_mul(lambda, aux1, aux2);    /* lambda = (3x + a) / 2y */

    mpz_mul(aux1, p_in->x, lambda);
    mpz_sub(nu, p_in->y, aux1);        /* nu = y - lambda * x */

    mpz_mul(aux1, lambda, lambda);
    mpz_submul_ui(aux1, p_in->x, 2);
    mpz_swap(p_out->x, aux1);          /* x' = lambda^2 - 2x */

    mpz_neg(p_out->y, nu);
    mpz_submul(p_out->y, p_out->x, lambda);   /* y' = -(lambda * x' + nu) */

    mpz_mod(p_out->x, p_out->x, ctx->m);
    mpz_mod(p_out->y, p_out->y, ctx->m);
clean:
#if !FAT_OBJECTS
    mpz_clears(aux1, aux2, lambda, nu, NULL);
#endif

    return !err;
}

int elliptic_mul(elliptic_point *p_out,
        elliptic_point *p_in,
        unsigned int times,
        elliptic_ctx *ctx)
{
    int err;
    elliptic_point point;
    elliptic_point_init_set(&point, p_in);

    while(times % 2 == 0)
    {
        err = elliptic_double(&point, &point, ctx);
        times >>= 1;

        if(err != 0) goto found;
    }

    elliptic_point_set(p_out, &point);

    err = elliptic_double(&point, &point, ctx);
    times >>= 1;

    if(err != 0) goto found;

    while(times)
    {
        if(times & 1)
        {
            err = elliptic_sum(p_out, p_out,
                    &point, ctx);
            if(err != 0) goto clean;
        }

        err = elliptic_double(&point, &point, ctx);
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
