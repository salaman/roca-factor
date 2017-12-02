#include <iostream>
#include <cstdio>
#include <chrono>
#include <vector>

#include <fplll.h>
#include <fmpq.h>
#include <fmpz.h>
#include <fmpz_poly.h>
#include <omp.h>

#include "parameters.h"
#include "discrete_log.h"

void print_fmpz_t(const char* fmt, const fmpz_t t)
{
    void (* freefunc)(void*, size_t);
    mp_get_memory_functions(nullptr, nullptr, &freefunc);

    char* str = fmpz_get_str(nullptr, 10, t);
    printf(fmt, str);
    freefunc(str, strlen(str) + 1);
}

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
void factor(const fmpz_t N, const fmpz_t gen, const fmpz_t Mp, const fmpz_t ord, unsigned int m, unsigned int t)
{
    fmpz_t cp, lowerBound, upperBound, cpPlusOrd;
    fmpz_init(cp);
    fmpz_init(lowerBound);
    fmpz_init(upperBound);
    fmpz_init(cpPlusOrd);

    // c' = log(N, 65537) mod M'
    discrete_log(cp, N, Mp, gen, ord);

    // Coppersmith parameters beta, X
    // beta: upper bounds
    // X
    // (b, X) <- (0.5, 2 * N^b/M')
    float b = 0.5;
    fmpz_t X;
    fmpz_init(X);
    fmpz_sqrt(X, N);
    fmpz_mul_ui(X, X, 2ul);
    fmpz_cdiv_q(X, X, Mp);

    // compute the lower bound of the search
    // c' / 2
    fmpz_fdiv_q_ui(lowerBound, cp, 2ul);
    fmpz_add(cpPlusOrd, cp, ord);

    // compute the upper bound of the search
    // (c' + ord') / 2
    fmpz_fdiv_q_ui(upperBound, cpPlusOrd, 2ul);

    // compute the range (upper - lower bounds) for threading purposes
    fmpz_t boundRange;
    fmpz_init_set(boundRange, upperBound);
    fmpz_sub(boundRange, boundRange, lowerBound);
    fmpz_add_ui(boundRange, boundRange, 1);

    #pragma omp critical
    {
        printf("parameters:\n");
        printf("N: "); fmpz_print(N); printf("\n");
        printf("M': "); fmpz_print(Mp); printf("\n");
        printf("gen: "); fmpz_print(gen); printf("\n");
        printf("c': "); fmpz_print(cp); printf("\n");
        printf("ord': "); fmpz_print(ord); printf("\n");
        printf("X: "); fmpz_print(X); printf("\n");
        printf("bounds: ["); fmpz_print(lowerBound); printf(", "); fmpz_print(upperBound); printf("]\n");
    }

    // N^b
    fmpz_t Npowb;
    fmpz_init_set(Npowb, N);
    fmpz_sqrt(Npowb, Npowb);

    // (M'^-1 mod N)
    fmpz_t Mpinv;
    fmpz_init(Mpinv);
    fmpz_invmod(Mpinv, Mp, N);

    // order of f(x) polynomial
    const unsigned int dd = 1;

    // number of vectors in matrix for LLL
    const unsigned int nn = dd * m + t;

    // partial polynomial f(x) = x + X
    fmpz_poly_t x;
    fmpz_poly_init(x);
    fmpz_poly_set_coeff_fmpz(x, 1, X);

    printf("\nspawning threads:\n");

    bool found = false;

    #pragma omp parallel
    {
        // thread setup - divide range over number of threads
        int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();

        fmpz_t boundRange, lowerBoundThread, upperBoundThread;
        fmpz_init_set(boundRange, upperBound);
        fmpz_sub(boundRange, boundRange, lowerBound);
        fmpz_add_ui(boundRange, boundRange, 1);
        fmpz_init(lowerBoundThread);
        fmpz_init(upperBoundThread);

        fmpz_mul_ui(lowerBoundThread, boundRange, this_thread);
        fmpz_fdiv_q_ui(lowerBoundThread, lowerBoundThread, num_threads);
        fmpz_add(lowerBoundThread, lowerBoundThread, lowerBound);

        fmpz_mul_ui(upperBoundThread, boundRange, this_thread + 1);
        fmpz_fdiv_q_ui(upperBoundThread, upperBoundThread, num_threads);
        fmpz_add(upperBoundThread, upperBoundThread, lowerBound);

        #pragma omp critical
        {
            printf("[thread %d] starting: [", this_thread); fmpz_print(lowerBoundThread); printf(", "); fmpz_print(upperBoundThread); printf(")\n");
        }

        // thread local init
        fmpz_t coeff, genpowap, ap, Npowm, Xpowi;
        fmpz_init(coeff);
        fmpz_init(genpowap);
        fmpz_init(Npowm);
        fmpz_init(Xpowi);

        fmpz_t divisor;
        fmpz_init(divisor);

        fmpz_poly_t newpol, fxpowi, fxpowm;
        fmpz_poly_init2(newpol, t);
        fmpz_poly_init(fxpowi);
        fmpz_poly_init(fxpowm);

        // partial polynomial f(x) = x
        fmpz_poly_t f;
        fmpz_poly_init2(f, 2);
        fmpz_poly_set_coeff_ui(f, 1, 1ul);

        fmpz_poly_t fx;
        fmpz_poly_init(fx);

        fmpz_t remainder, kp, result;
        fmpz_init(remainder);
        fmpz_init(kp);
        fmpz_init(result);
        // end thread local init

        unsigned long itr = 0;
        for (fmpz_init_set(ap, lowerBoundThread); fmpz_cmp(ap, upperBoundThread) < 1; fmpz_add_ui(ap, ap, 1ul), ++itr) {
            //auto start = std::chrono::high_resolution_clock::now();

            //#pragma omp critical
            //{
            //    printf("[thread %d] ap: ", this_thread); fmpz_print(ap); printf("\n");
            //}

            // f(x) <- x + (M'^-1 mod N) * (65537^a' mod M') mod N
            fmpz_powm(genpowap, gen, ap, Mp);

            fmpz_mul(coeff, Mpinv, genpowap);
            fmpz_mod(coeff, coeff, N);
            fmpz_poly_set_coeff_fmpz(f, 0, coeff);
            //print_fmpz_t("coeff: %s\n", coeff);

            // k' <- coppersmith(f(x), N, b, m, t, X)
            // adapted from https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/coppersmith.sage
            fmpz_poly_compose(fx, f, x);

            // compute polynomials
            fmpz_poly_struct ggt[nn];

            for (unsigned int i = 0; i < m; ++i) {
                for (unsigned int j = 0; j < dd; ++j) {
                    fmpz_poly_struct* curr = ggt + (i * dd + j);
                    fmpz_poly_init2(curr, m);

                    fmpz_poly_pow(curr, x, j);

                    fmpz_pow_ui(Npowm, N, m - i);
                    fmpz_poly_scalar_mul_fmpz(curr, curr, Npowm);

                    fmpz_poly_pow(fxpowi, fx, i);
                    fmpz_poly_mul(curr, curr, fxpowi);
                }
            }

            for (unsigned int i = 0; i < t; ++i) {
                fmpz_poly_struct* curr = ggt + (m * dd + i);
                fmpz_poly_init2(curr, m);
                fmpz_poly_pow(curr, x, i);

                fmpz_poly_pow(fxpowm, fx, m);

                fmpz_poly_mul(curr, curr, fxpowm);
            }

            // construct lattice B
            ZZ_mat<mpz_t> BB;
            BB.gen_zero(nn, nn);

            for (int i = 0; i < nn; ++i) {
                for (int j = 0; j <= i; ++j) {
                    auto ptr = fmpz_poly_get_coeff_ptr(ggt + i, j);
                    if (ptr != nullptr) {
                        fmpz_get_mpz(BB(i, j).get_data(), ptr);
                    }
                }

                fmpz_poly_clear(ggt + i);
            }

            // LLL
            fplll::lll_reduction(BB);

            // transform shortest vector in polynomial
            fmpz_poly_zero(newpol);
            for (unsigned int i = 0; i < nn; ++i) {
                fmpz_pow_ui(Xpowi, X, i);

                fmpz_t tmp;
                fmpz_init_set_readonly(tmp, BB(0, i).get_data());
                fmpz_fdiv_q(divisor, tmp, Xpowi);
                fmpz_clear_readonly(tmp);

                fmpz_poly_set_coeff_fmpz(newpol, i, divisor);
            }

            // factor polynomial
            fmpz_poly_factor_t factor;
            fmpz_poly_factor_init(factor);
            fmpz_poly_factor_zassenhaus(factor, newpol);

            // test roots
            for (unsigned int i = 0; i < factor->num; ++i) {
                //printf("factor: ");
                //fmpz_poly_print(factor->p + i);
                //flint_printf(", exp: %wd\n", factor->exp + i);

                // only interested in factors of form ax + b = 0
                if ((factor->p + i)->length != 2) {
                    continue;
                }

                // root = -b/a
                fmpz_fdiv_qr(kp, remainder, fmpz_poly_get_coeff_ptr(factor->p + i, 0), fmpz_poly_get_coeff_ptr(factor->p + i, 1));
                fmpz_neg(kp, kp);

                // test if root is an integer
                if (fmpz_is_zero(remainder)) {
                    //fmpz_poly_evaluate_fmpz(result, f, kp);
                    //fmpz_set(kp, f(quotient).evaluate()._fmpz());

                    // check the candidate root
                    // p <- k' * M' + (65537^a' mod M')
                    fmpz_mul(kp, kp, Mp);
                    fmpz_add(kp, kp, genpowap);

                    // test p|N for solution
                    if (fmpz_divisible(N, kp)) {
                        #pragma omp atomic write
                        found = true;

                        printf("[thread %d] *** FOUND! p = ", this_thread); fmpz_print(kp); printf("\n");
                        FILE* resfile = fopen("result.txt", "w");
                        fmpz_fprint(resfile, kp);
                        fclose(resfile);

                        break;
                    }
                }
            }

            fmpz_poly_factor_clear(factor);

            #pragma omp flush(found)
            if (found) {
                break;
            }

            // print progress every 10000 attempts
            if (itr % 10000 == 0) {
                fmpz_t remaining;
                fmpz_init(remaining);
                fmpz_sub(remaining, upperBoundThread, ap);
                #pragma omp critical
                {
                    printf("[thread %d] progress: ", this_thread);
                    fmpz_print(ap); printf("/"); fmpz_print(upperBoundThread);
                    printf(", remaining: "); fmpz_print(remaining); printf("\n");
                }
                fmpz_clear(remaining);
                itr = 0;
            }

            //printf("[thread %d] time taken: %lli\n", this_thread, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start));
        }

        #pragma omp critical
        printf("[thread %d] quitting\n", this_thread);

        // thread local cleanup
        fmpz_poly_clear(f);
        fmpz_poly_clear(fx);
        fmpz_clear(result);
        fmpz_clear(remainder);
        fmpz_clear(kp);
        fmpz_clear(divisor);
        fmpz_poly_clear(newpol);
        fmpz_clear(coeff);
        fmpz_clear(genpowap);
        fmpz_clear(ap);
        fmpz_clear(lowerBoundThread);
        fmpz_clear(upperBoundThread);
        fmpz_clear(boundRange);
        fmpz_poly_clear(fxpowi);
        fmpz_clear(Npowm);
        fmpz_clear(Xpowi);
        fmpz_poly_clear(fxpowm);
        // end thread local cleanup
    }

    fmpz_clear(Mpinv);
    fmpz_clear(cp);
    fmpz_clear(lowerBound);
    fmpz_clear(upperBound);
    fmpz_clear(cpPlusOrd);
}

int main()
{
    unsigned int m, t;
    fmpz_t N, Mp, generator, ord;
    fmpz_init(N);
    fmpz_init(Mp);
    fmpz_init(ord);

    // base prime generator is 65537
    fmpz_init_set_ui(generator, 65537ul);

    // input modulo N = pq
    fmpz_set_str(N, "944e13208a280c37efc31c3114485e590192adbb8e11c87cad60cdef0037ce99278330d3f471a2538fa667802ed2a3c44a8b7dea826e888d0aa341fd664f7fa7", 16);

    // find optimal parameters (M', ord', m, t) for the given N
    parameters(Mp, ord, &m, &t, N);

    factor(N, generator, Mp, ord, m, t);

    fmpz_clear(N);
    fmpz_clear(Mp);
    fmpz_clear(ord);
    fmpz_clear(generator);

    return 0;
}