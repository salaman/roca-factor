#include "parameters.h"

#include <cstdio>

void order(fmpz_t rot, const fmpz_t generator, const fmpz_t modulus)
{
    fmpz_set_ui(rot, 1);

    fmpz_t pow;
    fmpz_init(pow);
    while (true) {
        fmpz_powm(pow, generator, rot, modulus);

        if (fmpz_is_one(pow)) {
            fmpz_clear(pow);
            return;
        }

        fmpz_add_ui(rot, rot, 1);
    }
}

void algorithm2(fmpz_t res, const fmpz_t M, const fmpz_factor_t M_decomposition, const fmpz_t generator, const fmpz_t ordp)
{
    fmpz_set(res, M);

    fmpz_t ord;
    fmpz_init(ord);
    for (unsigned int i = 0; i < M_decomposition->num; ++i) {
        fmpz* p = M_decomposition->p + i;
        order(ord, generator, p);

        if (!fmpz_divisible(ordp, ord)) {
            // eliminate this prime from M'
            fmpz_fdiv_q(res, res, p);
        }
    }

    fmpz_clear(ord);
}

void parameters(fmpz_t M, fmpz_t ord, unsigned int* m, unsigned int* t, const fmpz_t N)
{
    mp_bitcnt_t bits = fmpz_bits(N);

    if (bits >= 512 && bits <= 960) {
        fmpz_primorial(M, 167);
        fmpz_set_str(ord, "2454106387091158800", 10);
        *m = 5;
        *t = 6;
    } else if (bits >= 962 && bits <= 1952) {
        fmpz_primorial(M, 353);
        fmpz_set_str(ord, "18008752815808067874770317538232798996000", 10);
        *m = 4;
        *t = 5;
    } else if (bits >= 1984 && bits < 3072) {
        fmpz_primorial(M, 701);
        fmpz_set_str(ord, "6d9bc8b22f6a31bd9622dfb3d3b93df4a7f37f1443ad5880236fd6ee5c5e5a80", 16);
        *m = 6;
        *t = 7;
    } else if (bits >= 3072 && bits <= 3936) {
        fmpz_primorial(M, 701);
        fmpz_set_str(ord, "6d9bc8b22f6a31bd9622dfb3d3b93df4a7f37f1443ad5880236fd6ee5c5e5a80", 10);
        *m = 25;
        *t = 26;
    } else if (bits >= 3968 && bits <= 4096) {
        fmpz_primorial(M, 1427);
        fmpz_set_str(ord, "33c07a6215a59c179c46f668d09de25aa7a394ec6d453871b4342d55115c8b8d2086d445a0b8ee950ee16fa5dae8ffd231191d87b5d80", 10);
        *m = 7;
        *t = 8;
    } else {
        // unsupported bit length
        return;
    }

    fmpz_t generator, Mcandidate, pv, candidate, ordDivCandidate, Mnew, eliminate;
    fmpz_init(Mcandidate);
    fmpz_init(pv);
    fmpz_init(candidate);
    fmpz_init(ordDivCandidate);
    fmpz_init(Mnew);
    fmpz_init(eliminate);

    fmpz_init_set_ui(generator, 65537ul);

    while (true) {
        double bestReward = 0.0;
        fmpz_t pv;
        fmpz_init_set_ui(pv, 0);
        ulong ev = 0;
        fmpz_set(Mcandidate, M);

        fmpz_factor_t Mfactors;
        fmpz_factor_init(Mfactors);
        fmpz_factor(Mfactors, M);

        fmpz_factor_t factors;
        fmpz_factor_init(factors);
        fmpz_factor(factors, ord);

        for (unsigned int i = 0; i < factors->num; ++i) {
            for (unsigned int exp = 1; exp <= factors->exp[i]; ++exp) {
                // candidate prime factor p^e for elimination
                fmpz_pow_ui(candidate, factors->p + i, exp);

                // ord/p^e
                fmpz_fdiv_q(ordDivCandidate, ord, candidate);

                // find M'_new
                algorithm2(Mnew, M, Mfactors, generator, ordDivCandidate);

                // calculate reward
                double reward = static_cast<double>((fmpz_flog_ui(ord, 2)) - fmpz_flog_ui(ordDivCandidate, 2)) / (fmpz_flog_ui(M, 2) - fmpz_flog_ui(Mnew, 2));

                if (reward > bestReward) {
                    bestReward = reward;
                    fmpz_set(pv, factors->p + i);
                    ev = exp;
                    fmpz_set(Mcandidate, Mnew);
                }
            }
        }

        fmpz_factor_clear(factors);
        fmpz_factor_clear(Mfactors);

        if (fmpz_is_zero(pv)) {
            // no candidates found
            break;
        }

        if (fmpz_flog_ui(Mcandidate, 2) < fmpz_flog_ui(N, 2) / 4) {
            // done
            //printf("done eliminating\n");
            break;
        }

        fmpz_set(M, Mcandidate);

        fmpz_pow_ui(eliminate, pv, ev);
        fmpz_fdiv_q(ord, ord, eliminate);

        //printf("eliminating "); fmpz_print(pv); printf("^%lu (reward %f)\n", ev, bestReward);
    }

    fmpz_clear(eliminate);
    fmpz_clear(pv);
    fmpz_clear(candidate);
    fmpz_clear(ordDivCandidate);
    fmpz_clear(Mnew);
    fmpz_clear(Mcandidate);
    fmpz_clear(generator);
}
