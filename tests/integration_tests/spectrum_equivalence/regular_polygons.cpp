#include <stdexcept>
#include "src/common/physical_optimization/OptimizationList.h"
#include "src/group/Group.h"
#include "src/spin_algebra/Multiplicity.h"
#include "tests/integration_tests/spectrum_equivalence/equivalence_check.h"
#include "tests/tools/OptimizationListsGenerator.h"

using namespace common::physical_optimization;

void initialize_regular_polygon_system(model::ModelInput& model, double J_value) {
    auto number_of_spins = model.getMults().size();
    auto J = model.addSymbol("J", J_value);
    for (int i = 0; i < number_of_spins; ++i) {
        model.assignSymbolToIsotropicExchange(J, i % number_of_spins, (i+1) % number_of_spins);
    }
}

void initialize_same_g_factor_for_regular_polygon_system(
    model::ModelInput& model, 
    double g_value) {
    auto g = model.addSymbol("g", g_value);
    for (int i = 0; i < model.getMults().size(); ++i) {
        model.assignSymbolToGFactor(g, i);
    }
} 

group::Group generate_Dn_group_for_polygon(size_t size) {
    group::Permutation generator_rotation(size);
    for (int i = 0; i < size; ++i) {
        generator_rotation[i] = (i + 1) % size;
    }
    group::Permutation generator_reflection(size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
    }
    return group::Group({group::Group::Dihedral, size}, 
        {generator_rotation, generator_reflection});
}

group::Group generate_first_C2_group_for_polygon(size_t size) {
    group::Permutation generator_reflection(size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[size - 1 - i] = i;
    }
    return group::Group(group::Group::S2, 
        {generator_reflection});
}

// TODO: this group (in combine with previous one) 
// leads to an error in the case of size=6 and ITO
group::Group generate_second_C2_group_for_polygon(size_t size) {
    if (size % 2 == 1) {
        throw std::invalid_argument("Cannot construct second C2 group for odd polygon");
    }
    group::Permutation generator_reflection(size);
    for (int i = 0; i < size; ++i) {
        generator_reflection[i] = (i + size / 2) % size;
    }
    return group::Group(group::Group::S2, 
        {generator_reflection});
}

class three_center_regular_polygon : public SpectrumFinalEquivalenceTest {};

#define group_triangle_diff_one group::Group({group::Group::Dihedral, 3}, {{2, 0, 1}, {0, 2, 1}})
#define group_triangle_diff_two group::Group({group::Group::Dihedral, 3}, {{1, 2, 0}, {1, 0, 2}})

TEST_P(three_center_regular_polygon, NoGFactors) {
    std::vector<double> js = {10, 17.17, 33};

    for (auto Jfirst : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_regular_polygon_system(model, Jfirst);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

TEST_P(three_center_regular_polygon, SameGFactors) {
    std::vector<double> js = {10, 17.17, 33};
    std::vector<double> gs = {2.1, 3.2};

    for (auto J_value : js) {
        for (auto g_value : gs) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_regular_polygon_system(model, J_value);
            initialize_same_g_factor_for_regular_polygon_system(model, g_value);

            runner::Runner runner_simple(model);

            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }      
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    three_center_regular_polygon,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_polygon(3),
            generate_first_C2_group_for_polygon(3),
            group_triangle_diff_one
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4}),
        std::vector<spin_algebra::Multiplicity>({5, 5, 5})
    )
),
spectrum_final_equivalence_test_name_generator
);

class four_center_regular_polygon : public SpectrumFinalEquivalenceTest {};

#define group_square_diff_one group::Group({group::Group::Dihedral, 4}, {{{1, 2, 3, 0}, {1, 0, 3, 2}}})
#define group_square_diff_two group::Group({group::Group::Dihedral, 4}, {{{3, 0, 1, 2}, {0, 3, 2, 1}}})

TEST_P(four_center_regular_polygon, NoGFactors) {
    std::vector<double> js = {10, -17.17, 33};

    for (auto J_value : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_regular_polygon_system(model, J_value);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

TEST_P(four_center_regular_polygon, SameGFactors) {
    std::vector<double> js = {10, -17.17, 33};
    std::vector<double> gs = {2.1, 3.2};

    for (auto J_value : js) {
        for (auto g_value : gs) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_regular_polygon_system(model, J_value);
            initialize_same_g_factor_for_regular_polygon_system(model, g_value);

            runner::Runner runner_simple(model);

            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }             
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    four_center_regular_polygon,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_polygon(4), 
            generate_first_C2_group_for_polygon(4),
            generate_second_C2_group_for_polygon(4),
            group_square_diff_one, 
            group_square_diff_two
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

class five_center_regular_polygon : public SpectrumFinalEquivalenceTest {};

TEST_P(five_center_regular_polygon, NoGFactors) {
    std::vector<double> js = {10, 17.17, 33};

    for (auto Jfirst : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_regular_polygon_system(model, Jfirst);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    five_center_regular_polygon,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_polygon(5),
            generate_first_C2_group_for_polygon(5)
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2}), 
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3})
    )
),
spectrum_final_equivalence_test_name_generator
);

class six_center_regular_polygon : public SpectrumFinalEquivalenceTest {};

TEST_P(six_center_regular_polygon, NoGFactors) {
    std::vector<double> js = {10, -17.17, 33};

    for (auto Jfirst : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_regular_polygon_system(model, Jfirst);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner, true);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    six_center_regular_polygon,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_polygon(6),
            generate_first_C2_group_for_polygon(6),
            // TODO: generate_second_C2_group_for_polygon(6) (in combine with previous one) 
            // leads to an error in the case of ITO
            group::Group(group::Group::S2, {{2, 1, 0, 5, 4, 3}})
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2}), 
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3, 3})
    )
),
spectrum_final_equivalence_test_name_generator
);

class seven_center_regular_polygon : public SpectrumFinalEquivalenceTest {};

TEST_P(seven_center_regular_polygon, NoGFactors) {
    std::vector<double> js = {10, 17.17, 33};

    for (auto Jfirst : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_regular_polygon_system(model, Jfirst);

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    seven_center_regular_polygon,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({
            generate_Dn_group_for_polygon(7),
            generate_first_C2_group_for_polygon(7)
        })
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2, 2})
    )
),
spectrum_final_equivalence_test_name_generator
);

class eight_center_regular_polygon : public SpectrumFinalEquivalenceTest {};

TEST_P(eight_center_regular_polygon, NoGFactors) {
    std::vector<double> js = {10, -17.17, 33};

    for (auto Jfirst : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_regular_polygon_system(model, Jfirst);

        OptimizationList reference(OptimizationList::LEX);
        reference
        .TzSort()
        .EliminatePositiveProjections()
        .Symmetrize(generate_Dn_group_for_polygon(8))
        .NonAbelianSimplify();

        runner::Runner runner_reference(model, reference);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_reference, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    eight_center_regular_polygon,
    ::testing::Combine(
    ::testing::Values(
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_first_C2_group_for_polygon(8)),
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_second_C2_group_for_polygon(8)),
        OptimizationList(OptimizationList::ITO)
            .TzSort()
            .TSquaredSort()
            .EliminateNonMininalProjections()
            .Symmetrize(generate_first_C2_group_for_polygon(8))
            .Symmetrize(generate_second_C2_group_for_polygon(8))
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3, 3, 3, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4, 4, 4, 4, 4, 4})
    )
),
spectrum_final_equivalence_test_name_generator
);
