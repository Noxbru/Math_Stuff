#include <stdio.h>
#include <gmp.h>

#include "mont.h"

int main(int argc, const char *argv[])
{
    mont_ctx mctx;
    mpz_t in1, in2;
    mpz_t mod;
    mpz_t out;

    if(argc != 4)
        return 1;

    mpz_init_set_str(in1, argv[1], 0);
    mpz_init_set_str(in2, argv[2], 0);
    mpz_init_set_str(mod, argv[3], 0);
    mpz_init(out);

    mont_init(&mctx, mod);

    mont_transform(in1, in1, &mctx);
    mont_transform(in2, in2, &mctx);

    mont_mul(out, in1, in2, &mctx);
    mont_inv_transform(out, out, &mctx);
    gmp_printf("%Zd\n",out);

    return 0;
}
