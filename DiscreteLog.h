#pragma once

#include <fmpz.h>
#include <fmpz_factor.h>
#include <fmpzxx.h>
#include <map>

class DiscreteLog
{
private:
    /**
     * Generator used for power calculations, set to 65537
     */
    fmpz_t generator_;

    /**
     * Multiplicative order of generator
     */
    fmpz_t generatorOrder_;

    /**
     * Prime decomposition of generator order
     * <prime, exponent>
     */
    std::map<int, int> generatorOrderDecomposition_;

    unsigned int maxPrimeNum(const fmpz_t N);
    void order(fmpz_t rot, const fmpz_t modulus);

    /**
     * Solves CRT for moduli and remainders
     *
     * x ≡ remainder mod modulus
     *
     * @param rot return value
     * @param moduli list of moduli
     * @param remainders list of remainders
     */
    void chineseRemainder(fmpz_t rot, const std::vector<flint::fmpzxx>& moduli, const std::vector<flint::fmpzxx>& remainders);
public:
    DiscreteLog();

    /**
     * Simple discrete logarithm for generator = b (65537)
     *
     * logb(a) mod M <=> b^k = a mod M
     *
     * @param rot return value
     * @param element a
     * @param m modulus M
     */
    void dlog(fmpz_t rot, const fmpz_t element, const fmpz_t m);
    void algorithm2(fmpz_t res, const fmpz_t m, const fmpz_factor_t mDecomposition, const fmpz_t order);
    void parameters(fmpz_t Mp, fmpz_t ord, unsigned int* m, unsigned int* t, const fmpz_t N);
};
