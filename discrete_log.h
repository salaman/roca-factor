#pragma once

#include <fmpz.h>
#include <fmpz_factor.h>
#include <fmpzxx.h>
#include <map>

/**
 * Simple discrete logarithm for generator = b (65537).
 *
 * logb(a) mod M <=> b^k ≡ a mod M
 *
 * @param rot return value
 * @param element a
 * @param m modulus M
 * @param generator base generator
 * @param generatorOrder multiplicative order of generator in Z/Zm
 */
void discrete_log(fmpz_t rot, const fmpz_t element, const fmpz_t m, const fmpz_t generator, const fmpz_t generatorOrder);

/**
 * Solves CRT for moduli and remainders.
 *
 * x ≡ remainder mod modulus
 *
 * @param rot return value
 * @param moduli list of moduli
 * @param remainders list of remainders
 */
void chinese_remainder(fmpz_t rot, const std::vector<flint::fmpzxx>& moduli, const std::vector<flint::fmpzxx>& remainders);