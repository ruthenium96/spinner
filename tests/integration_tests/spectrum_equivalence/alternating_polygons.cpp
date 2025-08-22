#include "src/spin_algebra/Multiplicity.h"
#include "tests/integration_tests/spectrum_equivalence/equivalence_check.h"
#include "tests/tools/OptimizationListsGenerator.h"

using namespace common::physical_optimization;

void initialize_alternating_polygon_system(model::ModelInput& model, double J_first_value, double J_second_value) {
    auto number_of_spins = model.getMults().size();
    auto J_first = model.addSymbol("J1", J_first_value);
    auto J_second = model.addSymbol("J2", J_second_value);
    for (int i = 0; i < number_of_spins; i+=2) {
        model.assignSymbolToIsotropicExchange(J_first, i % number_of_spins, (i+1) % number_of_spins);
        model.assignSymbolToIsotropicExchange(J_second, (i+1) % number_of_spins, (i+2) % number_of_spins);
    }
}

void initialize_same_g_factor_for_alternating_polygon_system(
    model::ModelInput& model, 
    double g_value) {
    auto g = model.addSymbol("g", g_value);
    for (int i = 0; i < model.getMults().size(); ++i) {
        model.assignSymbolToGFactor(g, i);
    }
} 

group::Group generate_Dn_half_group_for_alternating_polygon(size_t size) {
    group::Permutation generator_rotation(size);
    for (int i = 0; i < size; ++i) {
        generator_rotation[i] = (i + 2) % size;
    }
    group::Permutation generator_reflection(size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
    }
    return group::Group({group::Group::Dihedral, size / 2},
        {generator_rotation, generator_reflection});
}

group::Group generate_first_C2_group_for_alternating_polygon(size_t size) {
    group::Permutation generator_reflection(size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
    }
    return group::Group(group::Group::S2, 
        {generator_reflection});
}

group::Group generate_second_C2_group_for_alternating_polygon(size_t size) {
    if (size % 4 != 0) {
        throw std::invalid_argument("Cannot construct second C2 group for non-4n polygon");
    }
    group::Permutation generator_reflection(size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[i] = (i + size / 2) % size;
    }
    return group::Group(group::Group::S2, 
        {generator_reflection});
}

class four_center_alternating_polygon : public SpectrumFinalEquivalenceTest {};

TEST_P(four_center_alternating_polygon, NoGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (auto [Jfirst, Jsecond] : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_alternating_polygon_system(model, Jfirst, Jsecond);

        runner::Runner runner_simple(model);
        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

TEST_P(four_center_alternating_polygon, SameGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};
    std::vector<double> gs = {1.9, 3.1};

    for (auto [Jfirst, Jsecond] : js) {
        for (auto g_value : gs) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_alternating_polygon_system(model, Jfirst, Jsecond);
            initialize_same_g_factor_for_alternating_polygon_system(model, g_value);

            runner::Runner runner_simple(model);
            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);                
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    four_center_alternating_polygon,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_first_C2_group_for_alternating_polygon(4),
            generate_second_C2_group_for_alternating_polygon(4)
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4, 4}),
        std::vector<spin_algebra::Multiplicity>({5, 5, 5, 5})
    )
),
spectrum_final_equivalence_test_name_generator
);


class six_center_alternating_polygon : public SpectrumFinalEquivalenceTest {};

TEST_P(six_center_alternating_polygon, NoGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (auto [Jfirst, Jsecond] : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_alternating_polygon_system(model, Jfirst, Jsecond);

        runner::Runner runner_simple(model);
        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    six_center_alternating_polygon,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_half_group_for_alternating_polygon(6),
            generate_first_C2_group_for_alternating_polygon(6)
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3, 3})
    )
),
spectrum_final_equivalence_test_name_generator
);

class eight_center_alternating_polygon : public SpectrumFinalEquivalenceTest {};

TEST_P(eight_center_alternating_polygon, NoGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (auto [Jfirst, Jsecond] : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_alternating_polygon_system(model, Jfirst, Jsecond);
        
        OptimizationList reference(OptimizationList::LEX);
        reference
        .TzSort()
        .EliminatePositiveProjections()
        .Symmetrize(generate_Dn_half_group_for_alternating_polygon(8))
        .NonAbelianSimplify();

        runner::Runner runner_reference(model, reference);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_reference, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    eight_center_alternating_polygon,
    ::testing::Combine(
    ::testing::Values(
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_first_C2_group_for_alternating_polygon(8)),
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_second_C2_group_for_alternating_polygon(8)),
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_first_C2_group_for_alternating_polygon(8))
            .Symmetrize(generate_second_C2_group_for_alternating_polygon(8))
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3, 3, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4, 4, 4, 4, 4, 4})
    )
),
spectrum_final_equivalence_test_name_generator
);
