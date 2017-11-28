#include "DiscreteLog.h"

#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <flint/fmpz_factor.h>

DiscreteLog::DiscreteLog()
{
    fmpz_init_set_ui(generator_, 65537ul);
    fmpz_init_set_ui(generatorOrder_, 1201200ul);

    // {2: 4, 3: 2, 13: 1, 5: 2, 7: 1}
    generatorOrderDecomposition_ = {
            {2,  4},
            {3,  1},
            {5,  2},
            {7,  1},
            {11, 1},
            {13, 1},
            //{17, 1},
            //{23, 1},
            //{29, 1},
            //{37, 1},
            //{41, 1},
            //{53, 1},
            //{83, 1},
    };
}

void DiscreteLog::dlog(fmpz_t rot, const fmpz_t element, const fmpz_t m)
{
    fmpz_t powmod;
    fmpz_init(powmod);

    fmpz_powm(powmod, element, generatorOrder_, m);
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

    for (const auto& p : generatorOrderDecomposition_) {
        fmpz_set_ui(prime, p.first);
        fmpz_pow_ui(primeToPower, prime, p.second);
        fmpz_fdiv_q(orderDivPrimePower, generatorOrder_, primeToPower);
        fmpz_powm(gDash, generator_, orderDivPrimePower, m);
        fmpz_powm(hDash, element, orderDivPrimePower, m);
        bool found = false;

        fmpz_t i, gDashToI;
        fmpz_init(gDashToI);

        for (fmpz_init_set_ui(i, 0); fmpz_cmp(i, primeToPower) < 0; fmpz_add_ui(i, i, 1)) {
            fmpz_powm(gDashToI, gDash, i, m);
            if (fmpz_cmp(gDashToI, hDash) == 0) {
                flint::fmpzxx r;
                flint::fmpzxx m;
                fmpz_set(r._fmpz(), i);
                fmpz_set(m._fmpz(), primeToPower);
                remainders.push_back(r);
                moduli.push_back(m);
                found = true;
                break;
            }
        }

        fmpz_clear(i);
        fmpz_clear(gDashToI);

        if (!found) {
            throw std::runtime_error("No modulo found");
        }
    }

    chineseRemainder(rot, moduli, remainders);

    fmpz_clear(primeToPower);
    fmpz_clear(orderDivPrimePower);
    fmpz_clear(gDash);
    fmpz_clear(hDash);
}

void DiscreteLog::chineseRemainder(fmpz_t rot, const std::vector<flint::fmpzxx>& moduli, const std::vector<flint::fmpzxx>& remainders)
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

unsigned int DiscreteLog::maxPrimeNum(const fmpz_t N)
{
    mp_bitcnt_t bits = fmpz_bits(N);

    if (bits >= 512 && bits <= 960) {
        return 39;
    } else if (bits >= 962 && bits <= 1952) {
        return 71;
    } else if (bits >= 1984 && bits <= 3936) {
        return 126;
    } else if (bits >= 3968 && bits <= 4096) {
        return 225;
    }

    throw std::runtime_error("N invalid bit size");
}

void DiscreteLog::algorithm2(fmpz_t res, const fmpz_t M, const fmpz_factor_t M_decomposition, const fmpz_t ordp)
{
    fmpz_set(res, M);

    fmpz_t ord;
    fmpz_init(ord);
    for (unsigned int i = 0; i < M_decomposition->num; ++i) {
        fmpz* p = M_decomposition->p + i;
        order(ord, p);

        if (!fmpz_divisible(ordp, ord)) {
            // eliminate this prime from M'
            fmpz_fdiv_q(res, res, p);
        }
    }

    fmpz_clear(ord);
}

void DiscreteLog::order(fmpz_t rot, const fmpz_t modulus)
{
    fmpz_set_ui(rot, 1);

    fmpz_t pow;
    fmpz_init(pow);
    while (true) {
        fmpz_powm(pow, generator_, rot, modulus);

        if (fmpz_is_one(pow)) {
            fmpz_clear(pow);
            return;
        }

        fmpz_add_ui(rot, rot, 1);
    }
}

void DiscreteLog::parameters(fmpz_t Mp, fmpz_t ord, unsigned int* m, unsigned int* t, const fmpz_t N)
{
    auto P = maxPrimeNum(N);

    fmpz_primorial(Mp, P);

    // TODO
}
