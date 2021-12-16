#include <components/spectrum/SpectrumBuilder.h>

#include <armadillo>

#include "common/runner/Runner.h"
#include "gtest/gtest.h"

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

size_t size_of_spectrum_without_degeneracy(const Spectrum& spectrum) {
    size_t accumulator = 0;
    for (const auto& subspectrum : spectrum.blocks) {
        accumulator += subspectrum.raw_data.size();
    }
    return accumulator;
}

size_t size_of_spectrum_with_degeneracy(const Spectrum& spectrum) {
    size_t accumulator = 0;
    for (const auto& subspectrum : spectrum.blocks) {
        accumulator += subspectrum.raw_data.size() * subspectrum.properties.degeneracy;
    }
    return accumulator;
}

void expect_spectrum_equivalence(const Spectrum& first, const Spectrum& second) {
    for (size_t i = 0; i < first.blocks.size(); ++i) {
        EXPECT_EQ(first.blocks[i].raw_data, second.blocks[i].raw_data);
    }
}


TEST(matrix_and_spectrum_bulders, size_consistence_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_22_333_4444_23456_tzsort) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.TzSort();

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_2222_3333_4444_S2_S2_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 1, 2);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 2, 3);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 3, 0);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_2222_3333_4444_tzsort_S2_S2_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.TzSort();
        runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 1, 2);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 2, 3);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 3, 0);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444_S3_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 1, 2);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 2, 0);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444_tzsort_S3_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.TzSort();
        runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 1, 2);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 2, 0);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444_S3_symmetrize_nonabeliansimplify) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner.NonAbelianSimplify();

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 1, 2);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 2, 0);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
    }
}

TEST(
    matrix_and_spectrum_bulders,
    size_consistence_222_333_444_tzsort_S3_symmetrize_nonabeliansimplify) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.TzSort();
        runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner.NonAbelianSimplify();

        double J_value = 10;
        auto J = runner.modifySymbols().addSymbol("J", J_value);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 0, 1);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 1, 2);
        runner.modifySymbols().assignSymbolToIsotropicExchange(J, 2, 0);

        runner.InitializeSSquared();

        runner.BuildMatrices();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_with_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_matrix_without_degeneracy(
                runner.getMatrix(common::QuantityEnum::S_total_squared)));

        runner.BuildSpectra();

        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(runner.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_NE(
            runner.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner.getSpectrum(common::QuantityEnum::S_total_squared)));
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
//    for (size_t i = 0; i < mults.size(); ++i) {
//        for (size_t j = 0; j < mults.size(); ++j) {
//            js(i, j) = NAN;
//        }
//    }
//
//    runner.AddIsotropicExchange(js);
//
//    runner.BuildMatrices();
//
//    EXPECT_EQ(runner.getIndexConverter().get_total_space_size(), size_of_matrix_with_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
//    EXPECT_NE(runner.getIndexConverter().get_total_space_size(), size_of_matrix_without_degeneracy(runner.getMatrix(QuantityEnum::Energy)));
//}

TEST(spectrum_builder_without_matrix, size_consistence_and_spectra_equivalence_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(
    spectrum_builder_without_matrix,
    size_consistence_and_spectra_equivalence_22_333_4444_23456_tzsort) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        runner_without_matrices.TzSort();
        runner_using_matrices.TzSort();

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(
    spectrum_builder_without_matrix,
    size_consistence_and_spectra_equivalence_2222_3333_4444_S2_S2_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        runner_without_matrices.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        runner_without_matrices.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});
        runner_using_matrices.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        runner_using_matrices.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            1,
            2);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            2,
            3);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            3,
            0);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            1,
            2);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            2,
            3);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            3,
            0);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(
    spectrum_builder_without_matrix,
    size_consistence_and_spectra_equivalence_2222_3333_4444_tzsort_S2_S2_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        runner_without_matrices.TzSort();
        runner_without_matrices.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        runner_without_matrices.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});
        runner_using_matrices.TzSort();
        runner_using_matrices.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        runner_using_matrices.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            1,
            2);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            2,
            3);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            3,
            0);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            1,
            2);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            2,
            3);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            3,
            0);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(
    spectrum_builder_without_matrix,
    size_consistence_and_spectra_equivalence_222_333_444_S3_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        runner_without_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner_using_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            1,
            2);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            2,
            0);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            1,
            2);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            2,
            0);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(
    spectrum_builder_without_matrix,
    size_consistence_and_spectra_equivalence_222_333_444_tzsort_S3_symmetrize) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        runner_without_matrices.TzSort();
        runner_without_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner_using_matrices.TzSort();
        runner_using_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            1,
            2);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            2,
            0);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            1,
            2);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            2,
            0);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(
    spectrum_builder_without_matrix,
    size_consistence_and_spectra_equivalence_222_333_444_S3_symmetrize_nonabeliansimplify) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        runner_without_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner_without_matrices.NonAbelianSimplify();
        runner_using_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner_using_matrices.NonAbelianSimplify();

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            1,
            2);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            2,
            0);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            1,
            2);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            2,
            0);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_NE(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_NE(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(
    spectrum_builder_without_matrix,
    size_consistence_and_spectra_equivalence_222_333_444_tzsort_S3_symmetrize_nonabeliansimplify) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_using_matrices(mults);

        runner_without_matrices.TzSort();
        runner_without_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner_without_matrices.NonAbelianSimplify();
        runner_using_matrices.TzSort();
        runner_using_matrices.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        runner_using_matrices.NonAbelianSimplify();

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            1,
            2);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            2,
            0);
        auto J_using_matrices = runner_using_matrices.modifySymbols().addSymbol("J", J_value);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            0,
            1);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            1,
            2);
        runner_using_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_using_matrices,
            2,
            0);

        runner_without_matrices.InitializeSSquared();
        runner_using_matrices.InitializeSSquared();

        runner_using_matrices.BuildMatrices();

        runner_without_matrices.BuildSpectra();
        runner_using_matrices.BuildSpectra();

        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_NE(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::Energy)));
        EXPECT_EQ(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_with_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));
        EXPECT_NE(
            runner_without_matrices.getIndexConverter().get_total_space_size(),
            size_of_spectrum_without_degeneracy(
                runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared)));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            runner_using_matrices.getSpectrum(common::QuantityEnum::Energy));
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            runner_using_matrices.getSpectrum(common::QuantityEnum::S_total_squared));
    }
}

TEST(spectrum_builder_apply_to_entity, spectra_equivalence_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_without_matrices(mults);
        runner::Runner runner_to_manually_call_apply_to_entity(mults);

        double J_value = 10;
        auto J_without_matrices = runner_without_matrices.modifySymbols().addSymbol("J", J_value);
        runner_without_matrices.modifySymbols().assignSymbolToIsotropicExchange(
            J_without_matrices,
            0,
            1);
        auto J_to_manually_call_apply_to_entity =
            runner_to_manually_call_apply_to_entity.modifySymbols().addSymbol("J", J_value);
        runner_to_manually_call_apply_to_entity.modifySymbols().assignSymbolToIsotropicExchange(
            J_to_manually_call_apply_to_entity,
            0,
            1);

        runner_without_matrices.InitializeSSquared();
        runner_to_manually_call_apply_to_entity.InitializeSSquared();

        runner_to_manually_call_apply_to_entity.BuildMatrices();

        runner_without_matrices.BuildSpectra();

        SpectrumBuilder spectrum_builder;
        Spectrum manually_energy_spectrum = spectrum_builder.apply_to_energy(
            runner_to_manually_call_apply_to_entity.getMatrix(common::QuantityEnum::Energy));
        Spectrum manually_s_squared_spectrum =
            spectrum_builder.apply_to_non_energy(runner_to_manually_call_apply_to_entity.getMatrix(
                common::QuantityEnum::S_total_squared));

        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::Energy),
            manually_energy_spectrum);
        expect_spectrum_equivalence(
            runner_without_matrices.getSpectrum(common::QuantityEnum::S_total_squared),
            manually_s_squared_spectrum);
    }
}

TEST(
    spectrum_builder_apply_to_entity,
    throw_size_inconsistent_non_energy_matrix_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner_tz_sorted(mults);
        runner::Runner runner_not_tz_sorted_to_manually_call_apply_to_entity(mults);

        runner_tz_sorted.TzSort();

        double J_value = 10;
        auto J_tz_sorted = runner_tz_sorted.modifySymbols().addSymbol("J", J_value);
        runner_tz_sorted.modifySymbols().assignSymbolToIsotropicExchange(J_tz_sorted, 0, 1);
        auto J_not_tz_sorted_to_manually_call_apply_to_entity =
            runner_not_tz_sorted_to_manually_call_apply_to_entity.modifySymbols().addSymbol(
                "J",
                J_value);
        runner_not_tz_sorted_to_manually_call_apply_to_entity.modifySymbols()
            .assignSymbolToIsotropicExchange(
                J_not_tz_sorted_to_manually_call_apply_to_entity,
                0,
                1);

        runner_tz_sorted.InitializeSSquared();
        runner_not_tz_sorted_to_manually_call_apply_to_entity.InitializeSSquared();

        runner_tz_sorted.BuildMatrices();
        runner_not_tz_sorted_to_manually_call_apply_to_entity.BuildMatrices();

        runner_tz_sorted.BuildSpectra();

        SpectrumBuilder spectrum_builder;
        Spectrum manually_energy_spectrum = spectrum_builder.apply_to_energy(
            runner_not_tz_sorted_to_manually_call_apply_to_entity.getMatrix(
                common::QuantityEnum::Energy));
        EXPECT_THROW(
            spectrum_builder.apply_to_non_energy(
                runner_tz_sorted.getMatrix(common::QuantityEnum::S_total_squared)),
            std::length_error);
    }
}