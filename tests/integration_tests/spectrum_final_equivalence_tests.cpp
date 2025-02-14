#include <cstddef>
#include <optional>
#include "gtest/gtest.h"
#include "magic_enum.hpp"
#include "src/common/Quantity.h"
#include "src/common/physical_optimization/OptimizationList.h"
#include "src/common/runner/Runner.h"

using QuantumValues = std::array<std::optional<double>, magic_enum::enum_count<common::QuantityEnum>()>;

std::vector<QuantumValues> construct_final_vector(runner::Runner& runner) {
    std::vector<QuantumValues> vector;

    auto factory = runner.getDataStructuresFactories();
    auto degeneracy_vector = factory.createVector();

    std::array<std::optional<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>, 
        magic_enum::enum_count<common::QuantityEnum>()> values_vectors;

    for (const auto& quantity_enum_ : magic_enum::enum_values<common::QuantityEnum>()) {
        if (runner.getSpectrum(quantity_enum_).has_value()) {
            values_vectors[magic_enum::enum_integer<common::QuantityEnum>(quantity_enum_)] 
                = factory.createVector();
        } else {
            values_vectors[magic_enum::enum_integer<common::QuantityEnum>(quantity_enum_)] 
                = std::nullopt;
        }
    }

    for (const auto& subspectrum : runner.getSpectrum(common::Energy)->get().blocks) {
        degeneracy_vector->add_identical_values(
            subspectrum.raw_data->size(),
            subspectrum.properties.degeneracy);
    }

    for (int j = 0; j < values_vectors.size(); ++j) {
        const auto& quantity_enum_ = magic_enum::enum_value<common::QuantityEnum>(j);
        if (values_vectors[j].has_value()) {
            for (const auto& subspectrum : runner.getSpectrum(quantity_enum_)->get().blocks) {
                values_vectors[j]->get()->concatenate_with(subspectrum.raw_data);
            }
        }
    }

    values_vectors[magic_enum::enum_integer<common::QuantityEnum>(common::Energy)]->get()->subtract_minimum();

    for (size_t i = 0; i < degeneracy_vector->size(); ++i) {
        QuantumValues quantum_values;
        for (size_t j = 0; j < values_vectors.size(); ++j) {
            if (values_vectors[j].has_value()) {
                quantum_values[j] = values_vectors[j].value()->at(i);
            }
        }
        for (size_t j = 0; j < (size_t)degeneracy_vector->at(i); ++j) {
            vector.push_back(quantum_values);
        }
    }
    std::stable_sort(vector.begin(), vector.end());
    return vector;
}

void expect_final_vectors_equivalence(runner::Runner& simple, runner::Runner& second) {

    auto first_vector = construct_final_vector(simple);
    auto second_vector = construct_final_vector(second);

    QuantumValues quantum_values_sum_first;
    QuantumValues quantum_values_sum_second;

    for (int j = 0; j < quantum_values_sum_first.size(); ++j) {
        if (j == magic_enum::enum_integer<common::QuantityEnum>(common::Energy)) {
            continue;
        }
        if (first_vector.at(0)[j].has_value() && second_vector.at(0)[j].has_value()) {
            quantum_values_sum_first[j] = 0;
            quantum_values_sum_second[j] = 0;
        }
    }

    double last_energy = INFINITY;

    for (size_t i = 0; i < first_vector.size(); ++i) {
        // TODO: epsilon
        EXPECT_NEAR(first_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value(), 
            second_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value(), 1e-9);

        // TODO: values may differ, but the sums of values for degenerate eigenvectors should be the same
        if (std::abs(last_energy - first_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value()) > 1e-5) {
            for (int j = 0; j < quantum_values_sum_first.size(); ++j) {
                if (quantum_values_sum_first[j].has_value() && quantum_values_sum_second[j].has_value()) {
                    EXPECT_NEAR(quantum_values_sum_first[j].value(), quantum_values_sum_second[j].value(), 1e-3);
                    quantum_values_sum_first[j] = 0;
                    quantum_values_sum_second[j] = 0;
                }
            }
        }
        for (int j = 0; j < quantum_values_sum_first.size(); ++j) {
            const auto& quantity_enum_ = magic_enum::enum_value<common::QuantityEnum>(j);
            if (quantum_values_sum_first[j].has_value() && quantum_values_sum_second[j].has_value()) {
                quantum_values_sum_first[j].value() += first_vector[i][quantity_enum_].value();
                quantum_values_sum_second[j].value() += second_vector[i][quantity_enum_].value();  
            }
        }
        last_energy = first_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_two_center_system(
    model::ModelInput& model,
    double J_value) {
    auto J = model.addSymbol("J1", J_value);
    model.assignSymbolToIsotropicExchange(J, 0, 1);
}

TEST(spectrum_final_equivalence, two_center) {
    std::vector<std::vector<spin_algebra::Multiplicity>> multss = {
        {2, 2},
        {3, 3},
        {4, 4},
        {5, 5},
        {6, 6},
        {7, 7},
        {8, 8}};
    std::vector<double> js = {-50, -100, -80, -66, -40.1, -20, -10, -1.01, 1, 2, 10, 20, 33, 46, 50, 75, 100};
    for (const auto& mults : multss) {
        for (auto J_value : js) {
            group::Group group(group::Group::S2, {{1, 0}});

            model::ModelInput model(mults);
            initialize_two_center_system(model, J_value);

            runner::Runner runner_simple(model);

            // TZ_SORTER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort();
                runner::Runner runner_tz_sorted(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(group);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(group);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().Symmetrize(group);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList \
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList
                    .TzSort()
                    .EliminatePositiveProjections();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(group)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList 
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(group);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_five_center_mirror_symmetry_exchange_chain(
    model::ModelInput& model,
    double external,
    double internal) {
    auto Jexternal = model.addSymbol("J1", external);
    model.assignSymbolToIsotropicExchange(Jexternal, 0, 1)
        .assignSymbolToIsotropicExchange(Jexternal, 3, 4);

    auto Jinternal = model.addSymbol("J2", internal);
    model.assignSymbolToIsotropicExchange(Jinternal, 1, 2)
        .assignSymbolToIsotropicExchange(Jinternal, 2, 3);
}

TEST(spectrum_final_equivalence, five_center_mirror_symmetry_chain) {
    std::vector<std::vector<spin_algebra::Multiplicity>> multss = {
        {2, 2, 2, 2, 2},
        {2, 2, 3, 2, 2},
        {2, 2, 4, 2, 2},
        {3, 3, 2, 3, 3},
        {3, 3, 3, 3, 3},
        {3, 3, 4, 3, 3}};
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (const auto& mults : multss) {
        for (auto [Jfirst, Jsecond] : js) {
            group::Group group(group::Group::S2, {{4, 3, 2, 1, 0}});

            model::ModelInput model(mults);
            initialize_five_center_mirror_symmetry_exchange_chain(model, Jfirst, Jsecond);

            runner::Runner runner_simple(model);

            // TZ_SORTER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort();
                runner::Runner runner_tz_sorted(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(group);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(group);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().Symmetrize(group);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().TSquaredSort().EliminatePositiveProjections().SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList 
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList
                    .TzSort()
                    .EliminatePositiveProjections();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(group)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList 
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(group);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_four_centers_exchange_rectangle(
    model::ModelInput& model,
    double first,
    double second) {
    auto Jfirst = model.addSymbol("J1", first);
    model.assignSymbolToIsotropicExchange(Jfirst, 0, 1)
        .assignSymbolToIsotropicExchange(Jfirst, 2, 3);

    auto Jsecond = model.addSymbol("J2", second);
    model.assignSymbolToIsotropicExchange(Jsecond, 1, 2)
        .assignSymbolToIsotropicExchange(Jsecond, 3, 0);
}

TEST(spectrum_final_equivalence, rectangle) {
    std::vector<std::vector<spin_algebra::Multiplicity>> multss = {
        {2, 2, 2, 2},
        {3, 3, 3, 3},
        {4, 4, 4, 4}};
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (const auto& mults : multss) {
        for (auto [Jfirst, Jsecond] : js) {
            group::Group first_direction(group::Group::S2, {{1, 0, 3, 2}});
            group::Group second_direction(group::Group::S2, {{3, 2, 1, 0}});

            model::ModelInput model(mults);
            initialize_four_centers_exchange_rectangle(model, Jfirst, Jsecond);

            runner::Runner runner_simple(model);

            // TZ_SORTER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort();
                runner::Runner runner_tz_sorted(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(first_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(second_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(first_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(second_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(first_direction).Symmetrize(second_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(first_direction).Symmetrize(second_direction);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(first_direction)
                    .Symmetrize(second_direction);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList 
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList
                    .TzSort()
                    .EliminatePositiveProjections();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(first_direction)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(second_direction)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList 
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(first_direction);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                common::physical_optimization::OptimizationList 
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(second_direction);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + SYMMETRIZER + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(first_direction)
                    .Symmetrize(second_direction)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList 
                    optimizationList(common::physical_optimization::OptimizationList::ITO);
                optimizationList.TzSort()
                   .EliminatePositiveProjections()
                   .Symmetrize(first_direction)
                   .Symmetrize(second_direction);

               runner::Runner runner(model, optimizationList);
               expect_final_vectors_equivalence(runner_simple, runner);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_three_centers_exchange_triangle(model::ModelInput& model, double first) {
    auto Jfirst = model.addSymbol("J1", first);
    model.assignSymbolToIsotropicExchange(Jfirst, 0, 1)
        .assignSymbolToIsotropicExchange(Jfirst, 1, 2)
        .assignSymbolToIsotropicExchange(Jfirst, 2, 0);
}

TEST(spectrum_final_equivalence, triangle) {
    std::vector<std::vector<spin_algebra::Multiplicity>> multss = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};
    std::vector<double> js = {10, 17.17, 33};

    for (const auto& mults : multss) {
        for (auto Jfirst : js) {
            group::Group triangle(group::Group::S3, {{2, 0, 1}, {0, 2, 1}});
            model::ModelInput model(mults);
            initialize_three_centers_exchange_triangle(model, Jfirst);

            runner::Runner runner_simple(model);

            // TZ_SORTER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(triangle);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                group::Group triangle_diff(group::Group::S3, {{1, 2, 0}, {2, 1, 0}});
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(triangle_diff);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(triangle);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().Symmetrize(triangle);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER + NON_ABELIAN_SIMPLIFIER
            {}  // TZ_SORTER + SYMMETRIZER + NON_ABELIAN_SIMPLIFIER
            {}  // SYMMETRIZER + TZ_SORTER + NON_ABELIAN_SIMPLIFIER
            {}  // SYMMETRIZER + NON_ABELIAN_SIMPLIFIER + TZ_SORTER
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .TSquaredSort()
                    .EliminatePositiveProjections()
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // ITO + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList
                    .TzSort()
                    .EliminatePositiveProjections();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
        }
    }
}