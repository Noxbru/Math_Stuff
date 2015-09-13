#ifndef MONT_H
#define MONT_H

#include <gmp.h>
#include <stddef.h>

typedef struct mont_ctx
{
    mpz_t m;    /* The original modulus */
    mpz_t mm;
    size_t r;    /* The new modulus */
    mpz_t rr;
} mont_ctx;

/* The four numbers of the context must fulfill the relation:
 * m*mm - r*rr = 1
 */

void mont_init(mont_ctx *ctx, mpz_srcptr mod);

static inline void mont_transform(mpz_ptr out, mpz_srcptr n, mont_ctx *ctx)
{
    mpz_mul_2exp(out, n, ctx->r);
    mpz_mod(out, out, ctx->m);
}

static inline void mont_inv_transform(mpz_ptr out, mpz_srcptr n, mont_ctx *ctx)
{
    mpz_mul(out, n, ctx->rr);
    mpz_mod(out, out, ctx->m);
}

void mont_mul(mpz_ptr out, mpz_srcptr in1, mpz_srcptr in2, mont_ctx *ctx);

void mont_pow_ui(mpz_ptr out, mpz_srcptr in, unsigned long int times, mont_ctx *ctx);

void mont_pow_mpz(mpz_ptr out, mpz_srcptr in, mpz_srcptr times, mont_ctx *ctx);

#endif /* end of include guard: MONT_H */
