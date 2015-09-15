#include <stdio.h>
#include <gmp.h>

#include "mont.h"

int main(int argc, const char *argv[])
{
    mont_ctx mctx;
    mpz_t in;
    mpz_t times;
    mpz_t mod;
    mpz_t out;

    if(argc != 4)
        return 1;

    mpz_init_set_str(in, argv[1], 0);
    mpz_init_set_str(times, argv[2], 0);
    mpz_init_set_str(mod, argv[3], 0);
    mpz_init(out);

    mont_init(&mctx, mod);

    mont_transform(in, in, &mctx);

    mont_pow_mpz(out, in, times, &mctx);
    mont_inv_transform(out, out, &mctx);
    gmp_printf("%Zd\n",out);

    return 0;
}
