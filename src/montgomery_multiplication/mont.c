#include <assert.h>
#include <gmp.h>

#include "mont.h"

void mont_init(mont_ctx *ctx, mpz_srcptr mod)
{
    mpz_t not_used;
    mpz_t aux;

    mpz_init_set(ctx->m, mod);
    mpz_inits(ctx->mm, aux, ctx->rr, not_used, NULL);

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
    mpz_t aux1, aux2, aux3;

    mpz_inits(aux1, aux2, aux3, NULL);

    mpz_mul(aux1, in1, in2);
    mpz_mul(aux2, aux1, ctx->mm);
    mpz_tdiv_r_2exp(aux2, aux2, ctx->r);
    mpz_mul(aux3, aux2, ctx->m);

    mpz_add(aux1, aux1, aux3);
    mpz_tdiv_q_2exp(out, aux1, ctx->r);

    if(mpz_cmp(out, ctx->m) > 0)
        mpz_sub(out, out, ctx->m);
}
