#include <cstring>
#include <fmpz.h>

#include "parameters.h"
#include "factor.h"

int main(int argc, char** argv)
{
    if (argc < 2) {
        return 1;
    }

    unsigned int m, t;
    fmpz_t N, Mp, generator, ord;
    fmpz_init(N);
    fmpz_init(Mp);
    fmpz_init(ord);

    // base prime generator is 65537
    fmpz_init_set_ui(generator, 65537ul);

    // input modulo N = pq
    if (strncmp(argv[1], "0x", 2) == 0) {
        fmpz_set_str(N, argv[1] + 2, 16);
    } else {
        fmpz_set_str(N, argv[1], 10);
    }

    // find optimal parameters (M', ord', m, t) for the given N
    parameters(Mp, ord, &m, &t, N);

    factor(N, generator, Mp, ord, m, t);

    fmpz_clear(N);
    fmpz_clear(Mp);
    fmpz_clear(ord);
    fmpz_clear(generator);

    return 0;
}