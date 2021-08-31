#ifndef JULY_BASE_MATRICES_H
#define JULY_BASE_MATRICES_H

#include "old_code/indexes/blm.h"
#include <iostream>
#include <armadillo>

struct Matrices {
public:
    const Indexes & blm;

    // That is the strange way to use similar code for building matrices:
    virtual void add_to_hamiltonian(unsigned int i, unsigned int j, double v) = 0;
    virtual void add_to_s_squared(unsigned int i, unsigned int j, double v) = 0;
    typedef void (Matrices::*Add_To_Matrix)(unsigned int i, unsigned int j, double v);

    // Construct all matrices, do eigen-decomposition.
    virtual void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                                    arma::vec & degeneracy) = 0;

    explicit Matrices(const Indexes & blm): blm(blm) {};

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

    // "factor" is the factor before \sum_{i<j}(S_i; S_j)
    // for HDvV factor should be -2J
    // for S^2 factor should be 2:
    // S_{T}^2 = (\sum_{i}S_i)^2 = \sum_{i}S_i^2 + 2 * \sum_{i<j}(S_i; S_j)
    void scalar_product_total(Add_To_Matrix add_to_matrix, unsigned int repr, const double factor,
                              const unsigned long sym, const int a, const int b) {
        for (const Decomposition & d : blm.sym_to_block(sym)) {
            int na = blm.block_to_nzi(d.index, a);
            int nb = blm.block_to_nzi(d.index, b);

            // From block to sym and back again
            double decomposition_factor = d.coeff * blm.block_to_sym(d.index)[repr].coeff;

            // (Sa, Sb) = Sax*Sbx + Say*Sby + Saz*Sbz = 0.5 * (Sa+Sb- + Sa-Sb+) + Saz*Sbz

            // Saz Sbz
            (this->*add_to_matrix)(sym, sym, (na - blm.spins[a]) * (nb - blm.spins[b]) * factor * decomposition_factor);
            // Sa+ Sb-
            scalar_product_plus_minus(add_to_matrix, repr, factor, sym, d, a, b, na, nb);
            // Sa- Sb+
            scalar_product_plus_minus(add_to_matrix, repr, factor, sym, d, b, a, nb, na);
        }
    }
private:
    void scalar_product_plus_minus(Add_To_Matrix add_to_matrix, unsigned int repr, const double factor,
                                   const unsigned long sym, const Decomposition & d,
                                   const int a, const int b, const int na, const int nb) {
        if (na == blm.mults[a] - 1 || nb == 0) {
            return;
        }
        unsigned long new_block = blm.block_ladder(blm.block_ladder(d.index, a, +1), b, -1);

        if (repr < blm.block_to_sym(new_block).size()) {
            unsigned long new_sym = blm.block_to_sym(new_block)[repr].index;
            double decomposition_factor = d.coeff * blm.block_to_sym(new_block)[repr].coeff;

            // projection m = number n - spin S
            // so S(S+1)-m(m+1) = (2S-n)(n+1)
            // so S(S+1)-m(m-1) = n(2S+1-n)
            double factor_a = (2 * blm.spins[a] - na) * (na + 1);
            double factor_b = nb * (2 * blm.spins[b] + 1 - nb);

            (this->*add_to_matrix)(sym, new_sym, 0.5 * sqrt(factor_a * factor_b) * factor * decomposition_factor);
        }
    }
};

#endif