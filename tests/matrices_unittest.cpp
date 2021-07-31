#include "gtest/gtest.h"
#include <matrices/matrices.h>
#include <indexes/blm.h>

class Energy_S_Squared {
public:
    double energy;
    double s_squared;
    bool operator<(const Energy_S_Squared &other) const {
        return ((this->energy < other.energy) ||
        ((this->energy == other.energy) && (this->s_squared < other.s_squared)));
    }
};

// to unify all triples {energies, block_of_s_squared, degeneracy}
// returns 2-column matrix {energies, block_of_s_squared} sorted by energy
arma::dmat unify(arma::vec energy, arma::vec s_squared, arma::vec degeneracy) {
    EXPECT_EQ(energy.n_elem, s_squared.n_elem);
    EXPECT_EQ(energy.n_elem, degeneracy.n_elem);

    std::vector<Energy_S_Squared> vop;

    for (int i = 0; i < degeneracy.n_elem; ++i) {
        Energy_S_Squared ess; ess.energy = energy[i]; ess.s_squared = s_squared[i];
        for (int j = 0; j < degeneracy[i]; ++j) {
            vop.push_back(ess);
        }
    }
    std::sort(vop.begin(), vop.end());
    arma::dmat result(2, vop.size());
    for (int i = 0; i < vop.size(); ++i) {
        result(0, i) = vop[i].energy;
        result(1, i) = vop[i].s_squared;
    }
    return std::move(result);
}

TEST(Full_correctness, 22) {

    arma::vec eigval;
    arma::vec s_squared;
    arma::vec degeneracy;
    std::vector<int> mults = {2, 2};
    Indexes_P1 blm(mults);

    for (int i = 0; i < 100; ++i) {
        std::uniform_real_distribution<double> dist(0, 200.0);

        std::mt19937 rng;
        rng.seed(std::random_device{}());
        double J = dist(rng);


        arma::dmat js_arma(mults.size(), mults.size());
        js_arma(0, 1) = J; js_arma(1, 0) = J;

        Full_Matrices fmj(blm, js_arma);
        fmj.eigendecomposition(eigval, s_squared, degeneracy);
        EXPECT_EQ(2*J, eigval[3] - eigval[2]);
    }
}

TEST(Full_and_Sum_S_Square_equivalence, 2222) {

    arma::vec eigval;
    arma::vec s_squared;
    arma::vec degeneracy;
    std::vector<int> mults = {2, 2, 2, 2};
    Indexes_P1 blm(mults);
    blm.construct_branching_diagram();

    double J = 10;

    arma::dmat js_arma(mults.size(), mults.size());
    js_arma(0, 1) = J; js_arma(1, 0) = J;
    js_arma(1, 2) = J; js_arma(2, 1) = J;
    js_arma(2, 3) = J; js_arma(3, 2) = J;

    Full_Matrices fmj(blm, js_arma);
    fmj.eigendecomposition(eigval, s_squared, degeneracy);

    arma::dmat fmj_results = unify(eigval, s_squared, degeneracy);

    S2_Proj_Blocked_Matrices ssj(blm, js_arma);
    ssj.eigendecomposition(eigval, s_squared, degeneracy);

    arma::dmat ssj_results = unify(eigval, s_squared, degeneracy);

    arma::dmat diff = abs(fmj_results - ssj_results);

    EXPECT_TRUE(arma::approx_equal(ssj_results, fmj_results, "absdiff", 0.0001));
}