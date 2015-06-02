#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include <gmp.h>

typedef struct
{
    mpz_t A, B;
    mpz_t m;

#if FAT_OBJECTS
    mpz_t aux1, aux2;
    mpz_t lambda, nu;
#endif
} elliptic_ctx;

void static inline elliptic_init(elliptic_ctx *ctx)
{
    mpz_inits(ctx->A, ctx->B, ctx->m, NULL);
#if FAT_OBJECTS
    mpz_inits(ctx->aux1, ctx->aux2, ctx->lambda, ctx->nu, NULL);
#endif
}

void static inline elliptic_init_set(elliptic_ctx *ctx,
        mpz_ptr A, mpz_ptr B, mpz_ptr m)
{
    mpz_init_set(ctx->A, A);
    mpz_init_set(ctx->B, B);
    mpz_init_set(ctx->m, m);

#if FAT_OBJECTS
    mpz_inits(ctx->aux1, ctx->aux2, ctx->lambda, ctx->nu, NULL);
#endif
}

void static inline elliptic_clear(elliptic_ctx *ctx)
{
    mpz_clears(ctx->A, ctx->B, ctx->m, NULL);
#if FAT_OBJECTS
    mpz_clears(ctx->aux1, ctx->aux2, ctx->lambda, ctx->nu, NULL);
#endif
}

int elliptic_sum(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in1, mpz_ptr y_in1,
        mpz_ptr x_in2, mpz_ptr y_in2,
        elliptic_ctx *ctx);

int elliptic_double(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in, mpz_ptr y_in,
        elliptic_ctx *ctx);

int elliptic_mul(mpz_ptr x_out, mpz_ptr y_out,
        mpz_ptr x_in, mpz_ptr y_in,
        unsigned int times,
        elliptic_ctx *ctx);

#endif /* end of include guard: ELLIPTIC_H */
