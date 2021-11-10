#include "gtest/gtest.h"
#include "common/Runner.h"

#include <armadillo>

size_t size_of_matrix_without_degeneracy(const Matrix& matrix) {
    size_t accumulator = 0;
    for (const auto& submatrix : matrix.blocks) {
        accumulator += submatrix.raw_data.size();
    }
    return accumulator;
}

size_t size_of_matrix_with_degeneracy(const Matrix& matrix) {
    size_t accumulator = 0;
    for (const auto& submatrix : matrix.blocks) {
        accumulator += submatrix.raw_data.size() * submatrix.properties.degeneracy;
    }
    return accumulator;
}

size_t size_of_spectrum_without_degeneracy(const Matrix& matrix) {
    size_t accumulator = 0;
    for (const auto& submatrix : matrix.blocks) {
        accumulator += submatrix.raw_data.size();
    }
    return accumulator;
}

size_t size_of_spectrum_with_degeneracy(const Matrix& matrix) {
    size_t accumulator = 0;
    for (const auto& submatrix : matrix.blocks) {
        accumulator += submatrix.raw_data.size() * submatrix.properties.degeneracy;
    }
    return accumulator;
}

TEST(matrix_and_spectrum_bulders, throw_2222_inconsistent_symmetry) {
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
    js(0, 1) = 3*J; js(1, 0) = 3*J;
    js(1, 2) = J; js(2, 1) = J;
    js(2, 3) = 2*J; js(3, 2) = 2*J;
    js(3, 0) = J; js(0, 3) = J;

    runner.AddIsotropicExchange(js);

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(matrix_and_spectrum_bulders, size_consistence_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2},
                                                     {3, 3, 3},
                                                     {4, 4, 4, 4},
                                                     {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_22_333_4444_23456_tzsort) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2},
                                                     {3, 3, 3},
                                                     {4, 4, 4, 4},
                                                     {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.TzSort();

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_2222_3333_4444_S2_S2_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2},
                                                     {3, 3, 3, 3},
                                                     {4, 4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.Symmetrize(Group::S2, {{1, 0, 3, 2}});
        runner.Symmetrize(Group::S2, {{3, 2, 1, 0}});

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;
        js(1, 2) = J; js(2, 1) = J;
        js(2, 3) = J; js(3, 2) = J;
        js(3, 0) = J; js(0, 3) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_2222_3333_4444_tzsort_S2_S2_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2},
                                                     {3, 3, 3, 3},
                                                     {4, 4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
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
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;
        js(1, 2) = J; js(2, 1) = J;
        js(2, 3) = J; js(3, 2) = J;
        js(3, 0) = J; js(0, 3) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444_S3_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2},
                                                     {3, 3, 3},
                                                     {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.Symmetrize(Group::S3, {{1, 2, 0}, {0, 2, 1}});

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;
        js(1, 2) = J; js(2, 1) = J;
        js(2, 0) = J; js(0, 2) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444_tzsort_S3_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2},
                                                     {3, 3, 3},
                                                     {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.TzSort();
        runner.Symmetrize(Group::S3, {{1, 2, 0}, {0, 2, 1}});

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;
        js(1, 2) = J; js(2, 1) = J;
        js(2, 0) = J; js(0, 2) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444_S3_symmetrize_nonabeliansimplify) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2},
                                                     {3, 3, 3},
                                                     {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.Symmetrize(Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner.NonAbelianSimplify();

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;
        js(1, 2) = J; js(2, 1) = J;
        js(2, 0) = J; js(0, 2) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444_tzsort_S3_symmetrize_nonabeliansimplify) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2},
                                                     {3, 3, 3},
                                                     {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.TzSort();
        runner.Symmetrize(Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner.NonAbelianSimplify();

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;
        js(1, 2) = J; js(2, 1) = J;
        js(2, 0) = J; js(0, 2) = J;

        runner.AddIsotropicExchange(js);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
        EXPECT_EQ(runner.getTotalSpaceSize(), size_of_spectrum_with_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
        EXPECT_NE(runner.getTotalSpaceSize(), size_of_spectrum_without_degeneracy(runner.getMatrix(QuantityEnum::S_total_squared)));
    }
}

//TEST(matrix_and_spectrum_bulders, size_consistence_222222222_tzsort_S3_S3_symmetrize_nonabeliansimplify) {
//    std::vector<int> mults = {2, 2, 2,
//                              2, 2, 2,
//                              2, 2, 2};
//
//    runner::Runner runner(mults);
//
//    runner.TzSort();
//    runner.Symmetrize(Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
//    runner.NonAbelianSimplify();
//    runner.Symmetrize(Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
//    runner.NonAbelianSimplify();
//
//    arma::dmat js(mults.size(), mults.size());
//    for (int i = 0; i < mults.size(); ++i) {
//        for (int j = 0; j < mults.size(); ++j) {
//            js(i, j) = NAN;
//        }
//    }
//
//    runner.AddIsotropicExchange(js);
//
//    runner.BuildMatrices();
//
//    EXPECT_EQ(runner.getTotalSpaceSize(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
//    EXPECT_NE(runner.getTotalSpaceSize(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
//}