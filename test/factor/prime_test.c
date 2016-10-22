#include <stdio.h>
#include <stdlib.h>

#include "prime.h"

int main(int argc, const char *argv[])
{
    if(argc != 2)
        return 1;

    printf("%u\n",get_prime(atoi(argv[1])));

    return 0;
}
