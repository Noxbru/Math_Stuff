#include <stdio.h>
#include <gmp.h>

#include "factor.h"

int main(int argc, const char *argv[])
{
    mpz_t number, factor;

    if(argc != 2)
        return 1;

    mpz_init_set_str(number, argv[1], 0);
    mpz_init(factor);

    fermat(factor, number);
    gmp_printf("%Zd\n",factor);

    return 0;
}
