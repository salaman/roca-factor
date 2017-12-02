#pragma once

#include <fmpz.h>
#include <fmpz_factor.h>

/**
 * Order of generator in Z/Zm, ie. smallest non-zero integer such that generator^ord ≡ 1 mod M.
 *
 * @param rot return value
 * @param generator base
 * @param modulus M
 */
void order(fmpz_t rot, const fmpz_t generator, const fmpz_t modulus);

/**
 * The computation of the maximal divisor M' of the primorial M with ord_M'(65537)|ord' for a given ord' (divisor of ord_M(65537)).
 *
 * See "Algorithm 2" in paper.
 *
 * @param res return value
 * @param M modulus
 * @param M_decomposition prime decomposition of M
 * @param generator base
 * @param ordp order of generator in Z/ZM
 */
void algorithm2(fmpz_t res, const fmpz_t M, const fmpz_factor_t M_decomposition, const fmpz_t generator, const fmpz_t ordp);

/**
 * Returns optimized parameters to factor a given modulo N, based on its size.
 *
 * @param M return - optimized M'
 * @param ord return - order of generator in Z/ZM'
 * @param m return - optimal polynomial count for LLL
 * @param t return - optimal polynomial count for LLL
 * @param N modulo to factor
 */
void parameters(fmpz_t M, fmpz_t ord, unsigned int* m, unsigned int* t, const fmpz_t N);