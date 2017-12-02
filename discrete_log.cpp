#include "discrete_log.h"

#include <vector>
#include <stdexcept>
#include <cmath>
#include <flint/fmpz_factor.h>

void discrete_log(fmpz_t rot, const fmpz_t element, const fmpz_t m, const fmpz_t generator, const fmpz_t generatorOrder)
{
    fmpz_t powmod;
    fmpz_init(powmod);

    fmpz_powm(powmod, element, generatorOrder, m);
    if (fmpz_cmp_ui(powmod, 1) != 0) {
        throw std::runtime_error("Powmod not one");
    }

    std::vector<flint::fmpzxx> moduli;
    std::vector<flint::fmpzxx> remainders;

    fmpz_t prime, primeToPower, orderDivPrimePower, gDash, hDash;
    fmpz_init(prime);
    fmpz_init(primeToPower);
    fmpz_init(orderDivPrimePower);
    fmpz_init(gDash);
    fmpz_init(hDash);

    fmpz_factor_t generatorOrderFactors;
    fmpz_factor_init(generatorOrderFactors);
    fmpz_factor(generatorOrderFactors, generatorOrder);

    for (unsigned int i = 0; i < generatorOrderFactors->num; ++i) {
        fmpz_set(prime, generatorOrderFactors->p + i);
        fmpz_pow_ui(primeToPower, prime, generatorOrderFactors->exp[i]);
        fmpz_fdiv_q(orderDivPrimePower, generatorOrder, primeToPower);
        fmpz_powm(gDash, generator, orderDivPrimePower, m);
        fmpz_powm(hDash, element, orderDivPrimePower, m);
        bool found = false;

        fmpz_t j, gDashToI;
        fmpz_init(gDashToI);

        for (fmpz_init_set_ui(j, 0); fmpz_cmp(j, primeToPower) < 0; fmpz_add_ui(j, j, 1)) {
            fmpz_powm(gDashToI, gDash, j, m);
            if (fmpz_cmp(gDashToI, hDash) == 0) {
                flint::fmpzxx r;
                flint::fmpzxx m;
                fmpz_set(r._fmpz(), j);
                fmpz_set(m._fmpz(), primeToPower);
                remainders.push_back(r);
                moduli.push_back(m);
                found = true;
                break;
            }
        }

        fmpz_clear(j);
        fmpz_clear(gDashToI);

        if (!found) {
            throw std::runtime_error("No modulo found");
        }
    }

    chinese_remainder(rot, moduli, remainders);

    fmpz_clear(primeToPower);
    fmpz_clear(orderDivPrimePower);
    fmpz_clear(gDash);
    fmpz_clear(hDash);
}

void chinese_remainder(fmpz_t rot, const std::vector<flint::fmpzxx>& moduli, const std::vector<flint::fmpzxx>& remainders)
{
    if (moduli.size() != remainders.size()) {
        throw std::runtime_error("moduli and remainders are not of same size");
    }

    fmpz_t sum;
    fmpz_t product;
    fmpz_init_set_ui(sum, 0ul);
    fmpz_init_set_ui(product, 1ul);

    for (const auto& m : moduli) {
        fmpz_mul(product, product, m._fmpz());
    }

    fmpz_t p, pInvN, s;
    fmpz_init(p);
    fmpz_init(pInvN);
    fmpz_init(s);
    for (int i = 0; i < moduli.size(); ++i) {
        fmpz_fdiv_q(p, product, moduli.at(i)._fmpz());
        fmpz_invmod(pInvN, p, moduli.at(i)._fmpz());
        fmpz_mul(s, remainders.at(i)._fmpz(), pInvN);
        fmpz_mul(s, s, p);
        fmpz_add(sum, sum, s);
    }

    fmpz_mod(rot, sum, product);

    fmpz_clear(sum);
    fmpz_clear(product);
    fmpz_clear(p);
    fmpz_clear(pInvN);
    fmpz_clear(s);
}
