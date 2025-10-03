#include "gtest/gtest.h"
#include "src/group/Group.h"
#include "src/spin_algebra/Multiplicity.h"
#include "tests/integration_tests/spectrum_equivalence/equivalence_check.h"
#include "tests/tools/OptimizationListsGenerator.h"

using namespace common::physical_optimization;

void initialize_prism_system(model::ModelInput& model, double J_first_value, double J_second_value) {
    auto number_of_spins = model.getMults().size();

    auto J_first = model.addSymbol("J1", J_first_value);
    auto J_second = model.addSymbol("J2", J_second_value);

    for (int i = 0; i < number_of_spins / 2; ++i) {
        model.assignSymbolToIsotropicExchange(J_first, 
                                              i % (number_of_spins / 2), 
                                              (i+1) % (number_of_spins / 2));
        model.assignSymbolToIsotropicExchange(J_first, 
                                              i % (number_of_spins / 2) + number_of_spins / 2, 
                                              (i+1) % (number_of_spins / 2) + number_of_spins / 2);
        model.assignSymbolToIsotropicExchange(J_second, 
                                              i, 
                                              i + number_of_spins / 2);
    }
}

group::Group generate_Dn_group_for_prism(size_t size) {
    group::Permutation generator_rotation(2 * size);
    for (int i = 0; i < size; ++i) {
        generator_rotation[i] = (i + 1) % size;
        generator_rotation[i + size] = (i + 1) % size + size; 
    }
    group::Permutation generator_reflection(2 * size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
        generator_reflection[2 * size - 1 - i] = i + size;
    }
    return group::Group({group::Group::Dihedral, size}, 
        {generator_rotation, generator_reflection});
}

group::Group generate_C2_mirror_group_for_prism(size_t size) {
    group::Permutation generator_mirror(2 * size);
    for (int i = 0; i < size; ++i) {
        generator_mirror[i] = i + size;
        generator_mirror[i + size] = i;
    }
    return group::Group(group::Group::S2, 
        {generator_mirror});
}

group::Group generate_first_C2_group_for_prysm(size_t size) {
    group::Permutation generator_reflection(2 * size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
        generator_reflection[2 * size - 1 - i] = i + size;
    }
    return group::Group(group::Group::S2, 
        {generator_reflection});
}

class six_center_prism : public SpectrumFinalEquivalenceTest {};

TEST_P(six_center_prism, NoGFactors) {
    std::vector<std::pair<double, double>> js = {{10, -20}, {17.17, 15}, {33, -33}};

    for (auto jss : js) {
        auto Jfirst = jss.first;
        auto Jsecond = jss.second;
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_prism_system(model, Jfirst, Jsecond);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    six_center_prism,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_prism(3),
            generate_first_C2_group_for_prysm(3),
            generate_C2_mirror_group_for_prism(3)
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3, 3})
    )
),
spectrum_final_equivalence_test_name_generator
);