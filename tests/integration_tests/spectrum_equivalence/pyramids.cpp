#include "gtest/gtest.h"
#include "src/spin_algebra/Multiplicity.h"
#include "tests/integration_tests/spectrum_equivalence/equivalence_check.h"
#include "tests/tools/OptimizationListsGenerator.h"

using namespace common::physical_optimization;

void initialize_pyramid_system(model::ModelInput& model, double J_first_value, double J_second_value) {
    auto number_of_spins = model.getMults().size();
    auto J_first = model.addSymbol("J1", J_first_value);
    auto J_second = model.addSymbol("J2", J_second_value);
    for (int i = 0; i < number_of_spins - 1; ++i) {
        model.assignSymbolToIsotropicExchange(J_first, 
                                              i % (number_of_spins - 1), 
                                              (i+1) % (number_of_spins - 1));
        model.assignSymbolToIsotropicExchange(J_second, i, number_of_spins - 1);
    }
}

group::Group generate_Dn_group_for_pyramid(size_t size) {
    group::Permutation generator_rotation(size + 1);
    for (int i = 0; i < size; ++i) {
        generator_rotation[i] = (i + 1) % (size);
    }
    generator_rotation[size] = size;
    group::Permutation generator_reflection(size + 1);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
    }
    generator_reflection[size] = size;
    return group::Group({group::Group::Dihedral, size}, 
        {generator_rotation, generator_reflection});
}

group::Group generate_first_C2_group_for_pyramid(size_t size) {
    group::Permutation generator_reflection(size + 1);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
    }
    generator_reflection[size] = size;
    return group::Group(group::Group::S2, 
        {generator_reflection});
}

group::Group generate_second_C2_group_for_pyramid(size_t size) {
    if (size % 2 == 1) {
        throw std::invalid_argument("Cannot construct second C2 group for odd pyramid");
    }
    group::Permutation generator_reflection(size + 1);
    for (int i = 0; i < size; ++i) {
        generator_reflection[i] = (i + size / 2) % size;
    }
    generator_reflection[size] = size;
    return group::Group(group::Group::S2, 
        {generator_reflection});
}

class three_center_based_pyramid : public SpectrumFinalEquivalenceTest {};

TEST_P(three_center_based_pyramid, NoGFactors) {
    std::vector<std::vector<double>> jss = 
        {{10, -15}, {17.17, 11}, {33, 5}};

    for (const auto& js : jss) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_pyramid_system(model, js[0], js[1]);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    three_center_based_pyramid,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_pyramid(3),
            generate_first_C2_group_for_pyramid(3),
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

class four_center_based_pyramid : public SpectrumFinalEquivalenceTest {};

TEST_P(four_center_based_pyramid, NoGFactors) {
    std::vector<std::vector<double>> jss = 
        {{10, -15}, {17.17, 11}, {33, 5}};

    for (const auto& js : jss) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_pyramid_system(model, js[0], js[1]);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    four_center_based_pyramid,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_pyramid(4),
            generate_first_C2_group_for_pyramid(4),
            generate_second_C2_group_for_pyramid(4)
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4, 4, 4})
    )
),
spectrum_final_equivalence_test_name_generator
);

class six_center_based_pyramid : public SpectrumFinalEquivalenceTest {};

TEST_P(six_center_based_pyramid, NoGFactors) {
    std::vector<std::vector<double>> jss = 
        {{10, -15}, {17.17, 11}, {33, 5}};

    for (const auto& js : jss) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_pyramid_system(model, js[0], js[1]);

        OptimizationList reference(OptimizationList::LEX);
        reference
        .TzSort()
        .EliminatePositiveProjections()
        .Symmetrize(generate_Dn_group_for_pyramid(6))
        .NonAbelianSimplify();

        runner::Runner runner_reference(model, reference);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_reference, runner, true);
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    six_center_based_pyramid,
    ::testing::Combine(
    ::testing::Values(
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_first_C2_group_for_pyramid(6)),
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(group::Group(group::Group::S2, {{2, 1, 0, 5, 4, 3, 6}})),
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_first_C2_group_for_pyramid(6))
            // TODO: generate_second_C2_group_for_pyramid(6) (in combine with previous one) 
            // leads to an error in the case of ITO
            .Symmetrize(group::Group(group::Group::S2, {{2, 1, 0, 5, 4, 3, 6}}))
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3, 3, 2}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4, 4, 4, 4, 2})
    )
),
spectrum_final_equivalence_test_name_generator
);