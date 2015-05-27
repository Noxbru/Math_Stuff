#include <gmp.h>

#include "elliptic.h"

int elliptic_sum(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in1, mpz_ptr y_in1,
        mpz_ptr x_in2, mpz_ptr y_in2,
        elliptic_ctx ctx)
{
    int err;

    if(mpz_cmp(x_in1, x_in2) == 0)
    {
        if(mpz_cmp(y_in1, y_in2) == 0)
        {
            return elliptic_double(x_out, y_out, x_in1, y_in1, ctx);
        }
        else
        {
            // set the output point as infinity
            return 0;
        }
    }
    else
    {
        mpz_t aux1, aux2;
        mpz_t lambda, nu;

        mpz_inits(aux1, aux2, NULL);
        mpz_inits(lambda, nu, NULL);

        mpz_sub(aux1, x_in2, x_in1);
        err = mpz_invert(aux2, aux1, ctx.m);
        if(err == 0)
        {
            mpz_gcd(x_out, aux1, ctx.m);
            goto clean;
        }

        mpz_sub(aux1, y_in2, y_in1);
        mpz_mul(lambda, aux1, aux2);    /* lambda = (y2 - y1) / (x2 - x1) */

        mpz_mul(aux1, x_in1, lambda);
        mpz_sub(nu, y_in1, aux1);       /* nu = y1 - lambda * x1 */

        mpz_mul(aux1, lambda, lambda);
        mpz_sub(aux1, aux1, x_in1);
        mpz_sub(aux1, aux1, x_in2);
        mpz_swap(x_out, aux1);          /* x' = lambda^2 - x1 - x2 */

        mpz_neg(y_out, nu);
        mpz_submul(y_out, x_out, lambda);   /* y' = -(lambda * x' + nu) */

        mpz_mod(x_out, x_out, ctx.m);
        mpz_mod(y_out, y_out, ctx.m);

clean:
        mpz_clears(aux1, aux2, lambda, nu, NULL);
    }

    return !err;
}

int elliptic_double(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in, mpz_ptr y_in,
        elliptic_ctx ctx)
{
    int err;
    mpz_t aux1, aux2;
    mpz_t lambda, nu;

    if(mpz_cmp_ui(y_in, 0) == 0)
    {
        // set the output point as infinity
        return 0;
    }

    mpz_inits(aux1, aux2, NULL);
    mpz_inits(lambda, nu, NULL);

    mpz_add(aux1, y_in, y_in);
    err = mpz_invert(aux2, aux1, ctx.m);
    if(err == 0)
    {
        mpz_gcd(x_out, aux1, ctx.m);
        goto clean;
    }

    mpz_mul(aux1, x_in, x_in);
    mpz_mul_ui(aux1, aux1, 3);
    mpz_add(aux1, aux1, ctx.A);
    mpz_mul(lambda, aux1, aux2);    /* lambda = (3x + a) / 2y */

    mpz_mul(aux1, x_in, lambda);
    mpz_sub(nu, y_in, aux1);        /* nu = y - lambda * x */

    mpz_mul(aux1, lambda, lambda);
    mpz_submul_ui(aux1, x_in, 2);
    mpz_swap(x_out, aux1);          /* x' = lambda^2 - 2x */

    mpz_neg(y_out, nu);
    mpz_submul(y_out, x_out, lambda);   /* y' = -(lambda * x' + nu) */

    mpz_mod(x_out, x_out, ctx.m);
    mpz_mod(y_out, y_out, ctx.m);
clean:
    mpz_clears(aux1, aux2, lambda, nu, NULL);

    return !err;
}

int elliptic_mul(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in, mpz_ptr y_in,
        unsigned int times,
        elliptic_ctx ctx)
{
    int err;
    mpz_t aux1, aux2;
    mpz_inits(aux1, aux2, NULL);

    mpz_set(aux1, x_in);
    mpz_set(aux2, y_in);

    while(times % 2 == 0)
    {
        err = elliptic_double(aux1, aux2, aux1, aux2, ctx);
        times >>= 1;

        if(err != 0) goto found;
    }

    mpz_set(x_out, aux1);
    mpz_set(y_out, aux2);

    err = elliptic_double(aux1, aux2, aux1, aux2, ctx);
    times >>= 1;

    if(err != 0) goto found;

    while(times)
    {
        if(times & 1)
        {
            err = elliptic_sum(x_out, y_out, x_out, y_out,
                    aux1, aux2, ctx);
            if(err != 0) goto clean;
        }

        err = elliptic_double(aux1, aux2, aux1, aux2, ctx);
        times >>= 1;
        if(err != 0) goto found;
    }

clean:
    mpz_clears(aux1, aux2, NULL);
    return err;

found:
    mpz_set(x_out, aux1);
    goto clean;
}
