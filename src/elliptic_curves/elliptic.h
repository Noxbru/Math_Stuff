#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include <gmp.h>

/* Elliptic points >>> */
typedef enum
{
    POINT_AFFINE,
    POINT_PROJECTIVE,
    POINT_NULL
} point_type;

typedef struct
{
    mpz_t x, y, z;

    point_type type;
} elliptic_point;

void static inline elliptic_point_init(elliptic_point *point)
{
    mpz_inits(point->x, point->y, point->z, NULL);
    point->type = POINT_NULL;
}

void static inline elliptic_point_init_set(
        elliptic_point *point1,
        elliptic_point *point2)
{
    mpz_init_set(point1->x, point2->x);
    mpz_init_set(point1->y, point2->y);
    mpz_init_set(point1->z, point2->z);

    point1->type = point2->type;
}

void static inline elliptic_point_init_mpz_affine(
        elliptic_point *point,
        mpz_ptr x, mpz_ptr y)
{
    mpz_init_set(point->x, x);
    mpz_init_set(point->y, y);

    point->type = POINT_AFFINE;
}

void static inline elliptic_point_init_mpz_projective(
        elliptic_point *point,
        mpz_ptr x, mpz_ptr y, mpz_ptr z)
{
    mpz_init_set(point->x, x);
    mpz_init_set(point->y, y);
    mpz_init_set(point->z, z);

    point->type = POINT_PROJECTIVE;
}

void static inline elliptic_point_init_si_affine(
        elliptic_point *point,
        int x, int y)
{
    mpz_init_set_si(point->x, x);
    mpz_init_set_si(point->y, y);

    point->type = POINT_AFFINE;
}

void static inline elliptic_point_init_si_projective(
        elliptic_point *point,
        int x, int y, int z)
{
    mpz_init_set_si(point->x, x);
    mpz_init_set_si(point->y, y);
    mpz_init_set_si(point->z, z);

    point->type = POINT_PROJECTIVE;
}

void static inline elliptic_point_set(
        elliptic_point *point1,
        elliptic_point *point2)
{
    mpz_set(point1->x, point2->x);
    mpz_set(point1->y, point2->y);
    mpz_set(point1->z, point2->z);

    point1->type = point2->type;
}

void static inline elliptic_point_set_mpz_affine(
        elliptic_point *point,
        mpz_ptr x, mpz_ptr y)
{
    mpz_set(point->x, x);
    mpz_set(point->y, y);

    point->type = POINT_AFFINE;
}

void static inline elliptic_point_set_mpz_projective(
        elliptic_point *point,
        mpz_ptr x, mpz_ptr y, mpz_ptr z)
{
    mpz_set(point->x, x);
    mpz_set(point->y, y);
    mpz_set(point->z, z);

    point->type = POINT_PROJECTIVE;
}

void static inline elliptic_point_set_si_affine(
        elliptic_point *point,
        int x, int y)
{
    mpz_set_si(point->x, x);
    mpz_set_si(point->y, y);

    point->type = POINT_AFFINE;
}

void static inline elliptic_point_set_si_projective(
        elliptic_point *point,
        int x, int y, int z)
{
    mpz_set_si(point->x, x);
    mpz_set_si(point->y, y);
    mpz_set_si(point->z, z);

    point->type = POINT_PROJECTIVE;
}

void static inline elliptic_point_clear(elliptic_point *point)
{
    mpz_clears(point->x, point->y, point->z, NULL);
    point->type = POINT_NULL;
}

/* <<< */

/* Elliptic Curves >>> */
typedef enum
{
    ELLIPTIC_CURVE_WEIERSTRASS,
    ELLIPTIC_CURVE_NULL
} elliptic_curve_type;

typedef struct
{
    mpz_t A, B;
    mpz_t m;

#if FAT_OBJECTS
    mpz_t aux1, aux2;
    mpz_t lambda, nu;
#endif

    elliptic_curve_type type;
} elliptic_context;

/* Generic Elliptic Curves Functions >>> */
void static inline elliptic_curve_init(elliptic_context *ctx)
{
    mpz_inits(ctx->A, ctx->B, ctx->m, NULL);
#if FAT_OBJECTS
    mpz_inits(ctx->aux1, ctx->aux2, ctx->lambda, ctx->nu, NULL);
#endif
    ctx->type = ELLIPTIC_CURVE_NULL;
}

void static inline elliptic_curve_init_set(
        elliptic_context *ctx1,
        elliptic_context *ctx2)
{
    mpz_init_set(ctx1->A, ctx2->A);
    mpz_init_set(ctx1->B, ctx2->B);
    mpz_init_set(ctx1->m, ctx2->m);

#if FAT_OBJECTS
    mpz_inits(ctx1->aux1, ctx1->aux2, ctx1->lambda, ctx1->nu, NULL);
#endif

    ctx1->type = ctx2->type;
}

void static inline elliptic_curve_clear(elliptic_context *ctx)
{
    mpz_clears(ctx->A, ctx->B, ctx->m, NULL);
#if FAT_OBJECTS
    mpz_clears(ctx->aux1, ctx->aux2, ctx->lambda, ctx->nu, NULL);
#endif
    ctx->type = ELLIPTIC_CURVE_NULL;
}
/* <<< */



/* Weierstrass Elliptic Curves Functions >>> */
void static inline elliptic_curve_init_mpz_weierstrass(
        elliptic_context *ctx,
        mpz_ptr A, mpz_ptr B, mpz_ptr m)
{
    mpz_init_set(ctx->A, A);
    mpz_init_set(ctx->B, B);
    mpz_init_set(ctx->m, m);

#if FAT_OBJECTS
    mpz_inits(ctx->aux1, ctx->aux2, ctx->lambda, ctx->nu, NULL);
#endif

    ctx->type = ELLIPTIC_CURVE_WEIERSTRASS;
}

void static inline elliptic_curve_init_si_weierstrass(
        elliptic_context *ctx,
        int A, int B, int m)
{
    mpz_init_set_si(ctx->A, A);
    mpz_init_set_si(ctx->B, B);
    mpz_init_set_si(ctx->m, m);

#if FAT_OBJECTS
    mpz_inits(ctx->aux1, ctx->aux2, ctx->lambda, ctx->nu, NULL);
#endif

    ctx->type = ELLIPTIC_CURVE_WEIERSTRASS;
}

int elliptic_curve_sum_weierstrass_affine(
        elliptic_point *p_out,
        elliptic_point *p_in1,
        elliptic_point *p_in2,
        elliptic_context *ctx);

int elliptic_curve_double_weierstrass_affine(
        elliptic_point *p_out,
        elliptic_point *p_in,
        elliptic_context *ctx);

int elliptic_curve_mul_weierstrass_affine(
        elliptic_point *p_out,
        elliptic_point *p_in,
        unsigned int times,
        elliptic_context *ctx);

int elliptic_curve_sum_weierstrass_projective(
        elliptic_point *p_out,
        elliptic_point *p_in1,
        elliptic_point *p_in2,
        elliptic_context *ctx);

int elliptic_curve_double_weierstrass_projective(
        elliptic_point *p_out,
        elliptic_point *p_in,
        elliptic_context *ctx);

int elliptic_curve_mul_weierstrass_projective(
        elliptic_point *p_out,
        elliptic_point *p_in,
        unsigned int times,
        elliptic_context *ctx);
/* <<< */

/* Dispatching Elliptic Curves Functions >>> */
int static inline elliptic_curve_sum(
        elliptic_point *p_out,
        elliptic_point *p_in1,
        elliptic_point *p_in2,
        elliptic_context *ctx)
{
    if(ctx->type == ELLIPTIC_CURVE_WEIERSTRASS)
    {
        if(p_out->type == POINT_AFFINE)
            return elliptic_curve_sum_weierstrass_affine(p_out, p_in1, p_in2, ctx);
        /*else if(p_out->type == POINT_PROJECTIVE)*/
            /*return elliptic_curve_sum_weierstrass_projective(p_out, p_in1, p_in2, ctx);*/
    }
}

int static inline elliptic_curve_double(
        elliptic_point *p_out,
        elliptic_point *p_in,
        elliptic_context *ctx)
{
    if(ctx->type == ELLIPTIC_CURVE_WEIERSTRASS)
    {
        if(p_out->type == POINT_AFFINE)
            return elliptic_curve_double_weierstrass_affine(p_out, p_in, ctx);
        /*else if(p_out->type == POINT_PROJECTIVE)*/
            /*return elliptic_curve_double_weierstrass_projective(p_out, p_in, ctx);*/
    }
}

int static inline elliptic_curve_mul(
        elliptic_point *p_out,
        elliptic_point *p_in,
        unsigned int times,
        elliptic_context *ctx)
{
    if(ctx->type == ELLIPTIC_CURVE_WEIERSTRASS)
    {
        if(p_out->type == POINT_AFFINE)
            return elliptic_curve_mul_weierstrass_affine(p_out, p_in, times, ctx);
        /*else if(p_out->type == POINT_PROJECTIVE)*/
            /*return elliptic_curve_mul_weierstrass_projective(p_out, p_in, times, ctx);*/
    }
}
/* <<< */

#endif /* end of include guard: ELLIPTIC_H */
