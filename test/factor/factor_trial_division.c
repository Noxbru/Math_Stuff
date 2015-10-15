#include <stdio.h>
#include <gmp.h>

#include "factor.h"

int main(int argc, const char *argv[])
{
    mpz_t number;

    if(argc != 2)
        return 1;

    mpz_init_set_str(number, argv[1], 0);

    printf("%u\n", trial_division_ui(number));

    mpz_clear(number);

    return 0;
}
