#pragma once

#include <fmpz.h>

/**
 * Attempt factoring N
 *
 * @param N modulo N = pq
 * @param gen base generator
 * @param Mp optimized M'
 * @param ord multiplicative order of 65537 in Z/M'
 * @param m optimized number of f_i(x) polynomials (Coppersmith)
 * @param t optimized number of f_i+m(x) polynomials (Coppersmith)
 */
void factor(const fmpz_t N, const fmpz_t gen, const fmpz_t Mp, const fmpz_t ord, unsigned int m, unsigned int t);