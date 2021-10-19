#include "gtest/gtest.h"
#include "common/Runner.h"

#include <armadillo>

TEST(matrix_bulder, 22) {
    std::vector<int> mults = {2, 2};

    runner::Runner runner(mults);

    runner.TzSort();

    arma::dmat js(mults.size(), mults.size());
    for (int i = 0; i < mults.size(); ++i) {
        for (int j = 0; j < mults.size(); ++j) {
            js(i, j) = NAN;
        }
    }
    // 0 - 2
    // 3 - 1
    double J = 10;
    js(0, 1) = J; js(1, 0) = J;

    runner.AddIsotropicExchange(js);

    runner.BuildMatrix();

    std::cout << runner.getSpace() << std::endl;
}

TEST(matrix_bulder, 2222_S2_S2) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(Group::S2, {{3, 2, 1, 0}});


    arma::dmat js(mults.size(), mults.size());
    for (int i = 0; i < mults.size(); ++i) {
        for (int j = 0; j < mults.size(); ++j) {
            js(i, j) = NAN;
        }
    }
    // 0 - 2
    // 3 - 1
    double J = 10;
    js(0, 2) = J; js(2, 0) = J;
    js(1, 2) = J; js(2, 1) = J;
    js(1, 3) = J; js(3, 1) = J;
    js(3, 0) = J; js(0, 3) = J;

    runner.AddIsotropicExchange(js);

    runner.BuildMatrix();
}
