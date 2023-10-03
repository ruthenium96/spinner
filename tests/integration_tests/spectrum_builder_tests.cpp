#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

size_t size_of_matrix_without_degeneracy(const Matrix& matrix) {
    size_t accumulator = 0;
    for (const auto& submatrix : matrix.blocks) {
        accumulator += submatrix.raw_data->size();
    }
    return accumulator;
}

size_t size_of_matrix_with_degeneracy(const Matrix& matrix) {
    size_t accumulator = 0;
    for (const auto& submatrix : matrix.blocks) {
        accumulator += submatrix.raw_data->size() * submatrix.properties.degeneracy;
    }
    return accumulator;
}

size_t size_of_spectrum_without_degeneracy(const Spectrum& spectrum) {
    size_t accumulator = 0;
    for (const auto& subspectrum : spectrum.blocks) {
        accumulator += subspectrum.raw_data->size();
    }
    return accumulator;
}

size_t size_of_spectrum_with_degeneracy(const Spectrum& spectrum) {
    size_t accumulator = 0;
    for (const auto& subspectrum : spectrum.blocks) {
        accumulator += subspectrum.raw_data->size() * subspectrum.properties.degeneracy;
    }
    return accumulator;
}

void EXPECT_SIZE_CONSISTENCE_OF_MATRICES(runner::Runner& runner) {
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::Energy)));
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        size_of_matrix_with_degeneracy(runner.getMatrix(common::QuantityEnum::S_total_squared)));
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        size_of_matrix_without_degeneracy(runner.getMatrix(common::QuantityEnum::S_total_squared)));
}

void EXPECT_SIZE_CONSISTENCE_OF_SPECTRA(runner::Runner& runner) {
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

TEST(matrix_and_spectrum_bulders, size_consistence_22_333_4444_23456) {
    std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput model(mults);

        double J_value = 10;
        auto J = model.modifySymbolicWorker().addSymbol("J", J_value);
        model.modifySymbolicWorker().assignSymbolToIsotropicExchange(J, 0, 1);
        //
        {
            runner::Runner runner(model);

            runner.BuildSpectra();
            EXPECT_SIZE_CONSISTENCE_OF_MATRICES(runner);
            EXPECT_SIZE_CONSISTENCE_OF_SPECTRA(runner);
        }
        // TZ_SORT
        {
            common::physical_optimization::OptimizationList optimizationList;
            optimizationList.TzSort();
            runner::Runner runner(model, optimizationList);

            runner.BuildSpectra();
            EXPECT_SIZE_CONSISTENCE_OF_MATRICES(runner);
            EXPECT_SIZE_CONSISTENCE_OF_SPECTRA(runner);
        }
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_2222_3333_4444) {
    std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults = {
        {2, 2, 2, 2},
        {3, 3, 3, 3},
        {4, 4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput model(mults);
        double J_value = 10;
        auto J = model.modifySymbolicWorker().addSymbol("J", J_value);
        model.modifySymbolicWorker()
            .assignSymbolToIsotropicExchange(J, 0, 1)
            .assignSymbolToIsotropicExchange(J, 1, 2)
            .assignSymbolToIsotropicExchange(J, 2, 3)
            .assignSymbolToIsotropicExchange(J, 3, 0);

        // S2_S2_symmetrize
        {
            common::physical_optimization::OptimizationList optimizationList;
            optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
                .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

            runner::Runner runner(model, optimizationList);

            runner.BuildSpectra();
            EXPECT_SIZE_CONSISTENCE_OF_MATRICES(runner);
            EXPECT_SIZE_CONSISTENCE_OF_SPECTRA(runner);
        }
        // TZ_SORT + S2_S2_symmetrize
        {
            common::physical_optimization::OptimizationList optimizationList;
            optimizationList.TzSort()
                .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
                .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

            runner::Runner runner(model, optimizationList);

            runner.BuildSpectra();
            EXPECT_SIZE_CONSISTENCE_OF_MATRICES(runner);
            EXPECT_SIZE_CONSISTENCE_OF_SPECTRA(runner);
        }
    }
}

TEST(matrix_and_spectrum_bulders, size_consistence_222_333_444) {
    std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults = {
        {2, 2, 2},
        {3, 3, 3},
        {4, 4, 4}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput model(mults);
        double J_value = 10;
        auto J = model.modifySymbolicWorker().addSymbol("J", J_value);
        model.modifySymbolicWorker()
            .assignSymbolToIsotropicExchange(J, 0, 1)
            .assignSymbolToIsotropicExchange(J, 1, 2)
            .assignSymbolToIsotropicExchange(J, 2, 0);
        // S3_SYMMETRIZE
        {
            common::physical_optimization::OptimizationList optimizationList;
            optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});

            runner::Runner runner(model, optimizationList);

            runner.BuildSpectra();
            EXPECT_SIZE_CONSISTENCE_OF_MATRICES(runner);
            EXPECT_SIZE_CONSISTENCE_OF_SPECTRA(runner);
        }
        // TZ_SORT + S3_SYMMETRIZE
        {
            common::physical_optimization::OptimizationList optimizationList;
            optimizationList.TzSort().Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
            runner::Runner runner(model, optimizationList);

            runner.BuildSpectra();
            EXPECT_SIZE_CONSISTENCE_OF_MATRICES(runner);
            EXPECT_SIZE_CONSISTENCE_OF_SPECTRA(runner);
        }
        // S3_SYMMETRIZE + NON_ABELIAN_SIMPLIFY
        {}  // TZ_SORT + S3_SYMMETRIZE + NON_ABELIAN_SIMPLIFY
        {}
    }
}
