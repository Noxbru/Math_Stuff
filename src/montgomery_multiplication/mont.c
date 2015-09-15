#include <assert.h>
#include <gmp.h>

#include "mont.h"

void mont_init(mont_ctx *ctx, mpz_srcptr mod)
{
    mpz_t not_used;
    mpz_t aux;

    mpz_init_set(ctx->m, mod);
    mpz_inits(ctx->mm, aux, ctx->rr, not_used, NULL);

#if FAT_OBJECTS
    mpz_inits(ctx->aux1, ctx->aux2, ctx->aux3, NULL);
#endif

    ctx->r = mpz_sizeinbase(mod, 2);
    mpz_setbit(aux, mpz_sizeinbase(mod, 2));

    mpz_gcdext(not_used, ctx->mm, ctx->rr, ctx->m, aux);

    /* mpz_gcdext gives us numbers that satisfy the relation
     * m*mm + r*rr = 1 */
    mpz_neg(ctx->mm, ctx->mm);

    assert(mpz_cmp_ui(not_used, 1) == 0);
    mpz_clears(not_used, aux, NULL);
}

void mont_mul(mpz_ptr out, mpz_srcptr in1, mpz_srcptr in2, mont_ctx *ctx)
{
#if FAT_OBJECTS
    mpz_ptr aux1 = ctx->aux1;
    mpz_ptr aux2 = ctx->aux2;
    mpz_ptr aux3 = ctx->aux3;
#else
    mpz_t aux1, aux2, aux3;
    mpz_inits(aux1, aux2, aux3, NULL);
#endif

    mpz_mul(aux1, in1, in2);
    mpz_mul(aux2, aux1, ctx->mm);
    mpz_tdiv_r_2exp(aux2, aux2, ctx->r);
    mpz_mul(aux3, aux2, ctx->m);

    mpz_add(aux1, aux1, aux3);
    mpz_tdiv_q_2exp(out, aux1, ctx->r);

    if(mpz_cmp(out, ctx->m) > 0)
        mpz_sub(out, out, ctx->m);

#if !FAT_OBJECTS
    mpz_clears(aux1, aux2, aux3, NULL);
#endif
}

void mont_pow_ui(mpz_ptr out, mpz_srcptr in, unsigned long int times, mont_ctx *ctx)
{
    mpz_t aux;

    mpz_init_set_ui(aux, 1);
    mont_transform(aux, aux, ctx);
    mpz_set(out, in);

    while(times > 1)
    {
        if(times & 1)
            mont_mul(aux, aux, out, ctx);

        mont_mul(out, out, out, ctx);

        times >>=1;
    }

    mont_mul(out, out, aux, ctx);
}

void mont_pow_mpz(mpz_ptr out, mpz_srcptr in, mpz_srcptr times, mont_ctx *ctx)
{
    mpz_t aux, times_aux;

    if(mpz_fits_ulong_p(times))
    {
        mont_pow_ui(out, in, mpz_get_ui(times), ctx);
        return;
    }

    mpz_init_set_ui(aux, 1);
    mont_transform(aux, aux, ctx);
    mpz_init_set(times_aux, times);
    mpz_set(out, in);

    while(mpz_cmp_ui(times_aux, 1u) > 0)
    {
        if(mpz_tstbit(times_aux,0))
            mont_mul(aux, aux, out, ctx);

        mont_mul(out, out, out, ctx);

        mpz_tdiv_q_2exp(times_aux, times_aux, 1);
    }

    mont_mul(out, out, aux, ctx);
}
