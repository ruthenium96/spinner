#include "src/spin_algebra/Multiplicity.h"
#include "tests/integration_tests/spectrum_equivalence/equivalence_check.h"
#include "tests/tools/OptimizationListsGenerator.h"

using namespace common::physical_optimization;

void initialize_nine_centers_torus(model::ModelInput& model, double first) {
    auto J = model.addSymbol("J", first);
    model.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 3, 4)
        .assignSymbolToIsotropicExchange(J, 4, 5)
        .assignSymbolToIsotropicExchange(J, 6, 7)
        .assignSymbolToIsotropicExchange(J, 7, 8)
        .assignSymbolToIsotropicExchange(J, 0, 3)
        .assignSymbolToIsotropicExchange(J, 3, 6)
        .assignSymbolToIsotropicExchange(J, 1, 4)
        .assignSymbolToIsotropicExchange(J, 4, 7)
        .assignSymbolToIsotropicExchange(J, 2, 5)
        .assignSymbolToIsotropicExchange(J, 5, 8)
        .assignSymbolToIsotropicExchange(J, 2, 0)
        .assignSymbolToIsotropicExchange(J, 5, 3)
        .assignSymbolToIsotropicExchange(J, 8, 6)
        .assignSymbolToIsotropicExchange(J, 6, 0)
        .assignSymbolToIsotropicExchange(J, 7, 1)
        .assignSymbolToIsotropicExchange(J, 8, 2);
}

class torus : public SpectrumFinalEquivalenceTest {};

#define group_torus_hor group::Group({group::Group::Dihedral, 3}, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}})
#define group_torus_ver group::Group({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})

#define group_torus_hor_S2 group::Group(group::Group::S2, {{0, 1, 2, 6, 7, 8, 3, 4, 5}})
#define group_torus_ver_S2 group::Group(group::Group::S2, {{0, 2, 1, 3, 5, 4, 6, 8, 7}})

TEST_P(torus, NoGFactors) {
    std::vector<double> js = {10, 17.17, 33};

    for (auto j : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_nine_centers_torus(model, j);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner, true);
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    torus,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({group_torus_hor, group_torus_ver, group_torus_ver_S2, group_torus_hor_S2})
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2, 2, 2, 2})
    )
),
spectrum_final_equivalence_test_name_generator
);

