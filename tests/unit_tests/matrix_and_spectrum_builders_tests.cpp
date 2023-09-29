#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

// todo: do we really need this test?
TEST(
    spectrum_builder_apply_to_entity,
    throw_size_inconsistent_non_energy_matrix_22_333_4444_23456) {
    std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput model(mults);
        double J_value = 10;
        auto J_tz_sorted = model.modifySymbolicWorker().addSymbol("J", J_value);
        model.modifySymbolicWorker().assignSymbolToIsotropicExchange(J_tz_sorted, 0, 1);

        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.TzSort();

        runner::Runner runner_tz_sorted(model, optimizationList);
        runner::Runner runner_not_tz_sorted_to_manually_call_apply_to_entity(model);

        runner_tz_sorted.BuildMatrices();
        runner_not_tz_sorted_to_manually_call_apply_to_entity.BuildMatrices();

        runner_tz_sorted.BuildSpectra();

        auto [manually_energy_spectrum, unitary_transformation_matrices] =
            Spectrum::energy(runner_not_tz_sorted_to_manually_call_apply_to_entity.getMatrix(
                common::QuantityEnum::Energy));
        EXPECT_THROW(
            Spectrum::non_energy(
                runner_tz_sorted.getMatrix(common::QuantityEnum::S_total_squared),
                unitary_transformation_matrices),
            std::length_error);
    }
}