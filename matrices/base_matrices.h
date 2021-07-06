#ifndef JULY_BASE_MATRICES_H
#define JULY_BASE_MATRICES_H

#include "../indexes/blm.h"
#include <iostream>
#include <armadillo>

struct Matrices {
public:
    const Block_Lex_Map & blm;

    virtual void add_to_hamiltonian(unsigned int i, unsigned int j, double v) = 0;
    virtual void add_to_s_squared(unsigned int i, unsigned int j, double v) = 0;
    virtual void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                                    arma::vec & degeneracy) = 0;
    typedef void (Matrices::*Add_To_Matrix)(unsigned int i, unsigned int j, double v);

    explicit Matrices(const Block_Lex_Map & blm): blm(blm) {};

    void construct_mu_squared(const arma::vec & eigval, const arma::vec & s_squared,
                              const arma::vec & degeneracities) {
        double g = 2.0;
        arma::vec wep = s_squared % degeneracities;
        std::cout << "{" << std::endl;
        for (int t = 1; t < 301; ++t) {
            double mu_squared = g * g * sum(exp(- eigval / t) % wep);
            double statsum = sum(exp(- eigval / t) % degeneracities);
            std::cout << t << " : " << mu_squared / statsum << "," << std::endl;
        }
        std::cout << "}" << std::endl;
    }

    // build block of s_squared matrix from (start, start) to (end - 1, end - 1)
    void construct_block_of_s_squared(unsigned long start, unsigned long end, Add_To_Matrix s_squared_) {
        double diagonal_sum = 0;
        for (double s : blm.spins) {
            diagonal_sum += s * (s + 1);
        }
        #pragma omp parallel for
        for (unsigned long block = start; block < end; ++block) {
            add_to_s_squared(block, block, diagonal_sum);
            for (int a = 0; a < blm.v_size; ++a) {
                for (int b = a + 1; b < blm.v_size; ++b) {
                    if (a != b) {
                        scalar_product_total(s_squared_, 2, block, a, b);
                    }
                }
            }
        }
    }

    // "factor" is the factor before \sum_{i<j}(S_i; S_j)
    // for HDvV factor should be -2J
    // for S^2 factor should be 2:
    // S_{T}^2 = (\sum_{i}S_i)^2 = \sum_{i}S_i^2 + 2 * \sum_{i<j}(S_i; S_j)
    void scalar_product_total(Add_To_Matrix add_to_matrix, const double factor,
                              const unsigned long block, const int a, const int b) {
        int na = blm.block_to_nzi(block, a);
        int nb = blm.block_to_nzi(block, b);

        // Saz Sbz
        (this->*add_to_matrix)(block, block, (na - blm.spins[a]) * (nb - blm.spins[b]) * factor);
        // Sa+ Sb-
        scalar_product_plus_minus(add_to_matrix, factor, block, a, b, na, nb);
        // Sa- Sb+
        scalar_product_plus_minus(add_to_matrix, factor, block, b, a, nb, na);
    }
private:
    void scalar_product_plus_minus(Add_To_Matrix add_to_matrix, const double factor, const unsigned long block,
                                   const int a, const int b, const int na, const int nb) {
        if (na == blm.mults[a] - 1 || nb == 0) {
            return;
        }
        unsigned long new_block = blm.block_ladder(blm.block_ladder(block, a, +1), b, -1);
        // projection m = number n - spin S
        // so S(S+1)-m(m+1) = (2S-n)(n+1)
        // so S(S+1)-m(m-1) = n(2S+1-n)
        double factor_a = (2 * blm.spins[a] - na) * (na + 1);
        double factor_b = nb * (2 * blm.spins[b] + 1 - nb);
        (this->*add_to_matrix)(block, new_block, 0.5 * sqrt(factor_a * factor_b) * factor);
    }
};

#endif