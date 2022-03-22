#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

struct EnergyAndSSquared {
    bool operator<(const EnergyAndSSquared& rhs) const {
        if (energy < rhs.energy)
            return true;
        if (rhs.energy < energy)
            return false;
        return s_squared < rhs.s_squared;
    }
    bool operator>(const EnergyAndSSquared& rhs) const {
        return rhs < *this;
    }
    bool operator<=(const EnergyAndSSquared& rhs) const {
        return !(rhs < *this);
    }
    bool operator>=(const EnergyAndSSquared& rhs) const {
        return !(*this < rhs);
    }
    double energy;
    double s_squared;
};

std::vector<EnergyAndSSquared> construct_final_vector(const runner::Runner& runner) {
    std::vector<EnergyAndSSquared> vector;

    DenseVector energy_vector;
    DenseVector s_squared_vector;
    DenseVector degeneracy_vector;

    for (const auto& subspectrum : runner.getSpectrum(common::Energy).blocks) {
        energy_vector.concatenate_with(subspectrum.raw_data);
        degeneracy_vector.add_identical_values(
            subspectrum.raw_data.size(),
            subspectrum.properties.degeneracy);
    }
    energy_vector.subtract_minimum();

    for (const auto& subspectrum : runner.getSpectrum(common::S_total_squared).blocks) {
        s_squared_vector.concatenate_with(subspectrum.raw_data);
    }

    for (size_t i = 0; i < degeneracy_vector.size(); ++i) {
        EnergyAndSSquared energy_and_s_squared = {energy_vector(i), s_squared_vector(i)};
        for (size_t j = 0; j < (size_t)degeneracy_vector(i); ++j) {
            vector.push_back(energy_and_s_squared);
        }
    }
    std::stable_sort(vector.begin(), vector.end());
    return vector;
}

void expect_final_vectors_equivalence(const runner::Runner& simple, runner::Runner& second) {
    second.BuildSpectra();

    auto first_vector = construct_final_vector(simple);
    auto second_vector = construct_final_vector(second);

    double s_squared_sum_first = 0;
    double s_squared_sum_second = 0;
    double last_energy = INFINITY;

    for (size_t i = 0; i < first_vector.size(); ++i) {
        // TODO: epsilon
        EXPECT_NEAR(first_vector[i].energy, second_vector[i].energy, 1e-9);
        // TODO: now s2-squared values can be incorrect,
        // EXPECT_NEAR(first_vector[i].s_squared, second_vector[i].s_squared, 1e-9);

        // TODO: but sum of s_squared for degenerate eigenvectors is correct
        if (std::abs(last_energy - first_vector[i].energy) > 1e-5) {
            EXPECT_NEAR(s_squared_sum_first, s_squared_sum_second, 1e-3);
            s_squared_sum_first = 0;
            s_squared_sum_second = 0;
        }
        s_squared_sum_first += first_vector[i].s_squared;
        s_squared_sum_second += second_vector[i].s_squared;
        last_energy = first_vector[i].energy;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_four_centers_exchange_rectangle(model::Model& model, double first, double second) {
    auto Jfirst = model.getSymbols().addSymbol("J1", first);
    model.getSymbols()
        .assignSymbolToIsotropicExchange(Jfirst, 0, 1)
        .assignSymbolToIsotropicExchange(Jfirst, 2, 3);

    auto Jsecond = model.getSymbols().addSymbol("J2", second);
    model.getSymbols()
        .assignSymbolToIsotropicExchange(Jsecond, 1, 2)
        .assignSymbolToIsotropicExchange(Jsecond, 3, 0);

    model.InitializeSSquared();
}

TEST(spectrum_final_equivalence, rectangle) {
    std::vector<std::vector<int>> multss = {{2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (const auto& mults : multss) {
        for (auto [Jfirst, Jsecond] : js) {
            group::Group first_direction(group::Group::S2, {{1, 0, 3, 2}});
            group::Group second_direction(group::Group::S2, {{3, 2, 1, 0}});

            model::Model model(mults);
            initialize_four_centers_exchange_rectangle(model, Jfirst, Jsecond);

            runner::Runner runner_simple(model);

            runner_simple.BuildSpectra();

            // TZ_SORTER
            {
                runner::Runner runner_tz_sorted(model);
                runner_tz_sorted.TzSort();
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                runner::Runner runner_tz_sorted_eliminated(model);
                runner_tz_sorted_eliminated.TzSort();
                runner_tz_sorted_eliminated.EliminatePositiveProjections();
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted_eliminated);
            }
            // SYMMETRIZER
            {
                runner::Runner runner_first_direction(model);
                runner_first_direction.Symmetrize(first_direction);
                expect_final_vectors_equivalence(runner_simple, runner_first_direction);
            }
            {
                runner::Runner runner_second_direction(model);
                runner_second_direction.Symmetrize(second_direction);
                expect_final_vectors_equivalence(runner_simple, runner_second_direction);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                runner::Runner runner_first_direction_tz_sorted(model);
                runner_first_direction_tz_sorted.TzSort();
                runner_first_direction_tz_sorted.Symmetrize(first_direction);
                expect_final_vectors_equivalence(runner_simple, runner_first_direction_tz_sorted);
            }
            {
                runner::Runner runner_second_direction_tz_sorted(model);
                runner_second_direction_tz_sorted.TzSort();
                runner_second_direction_tz_sorted.Symmetrize(second_direction);
                expect_final_vectors_equivalence(runner_simple, runner_second_direction_tz_sorted);
            }
            // SYMMETRIZER + TZ_SORTER
            {
                runner::Runner runner_first_direction_tz_sorted(model);
                runner_first_direction_tz_sorted.Symmetrize(first_direction);
                runner_first_direction_tz_sorted.TzSort();
                expect_final_vectors_equivalence(runner_simple, runner_first_direction_tz_sorted);
            }
            {
                runner::Runner runner_second_direction_tz_sorted(model);
                runner_second_direction_tz_sorted.Symmetrize(second_direction);
                runner_second_direction_tz_sorted.TzSort();
                expect_final_vectors_equivalence(runner_simple, runner_second_direction_tz_sorted);
            }
            // SYMMETRIZER + SYMMETRIZER
            {
                runner::Runner runner_both_directions(model);
                runner_both_directions.Symmetrize(first_direction);
                runner_both_directions.Symmetrize(second_direction);
                expect_final_vectors_equivalence(runner_simple, runner_both_directions);
            }
            // TZ_SORTER + SYMMETRIZER + SYMMETRIZER
            {
                runner::Runner runner_both_directions_tz_sorted(model);
                runner_both_directions_tz_sorted.TzSort();
                runner_both_directions_tz_sorted.Symmetrize(first_direction);
                runner_both_directions_tz_sorted.Symmetrize(second_direction);
                expect_final_vectors_equivalence(runner_simple, runner_both_directions_tz_sorted);
            }
            // SYMMETRIZER + TZ_SORTER + SYMMETRIZER
            {
                runner::Runner runner_both_directions_tz_sorted(model);
                runner_both_directions_tz_sorted.Symmetrize(first_direction);
                runner_both_directions_tz_sorted.TzSort();
                runner_both_directions_tz_sorted.Symmetrize(second_direction);
                expect_final_vectors_equivalence(runner_simple, runner_both_directions_tz_sorted);
            }
            // SYMMETRIZER + SYMMETRIZER + TZ_SORTER
            {
                runner::Runner runner_both_directions_tz_sorted(model);
                runner_both_directions_tz_sorted.Symmetrize(first_direction);
                runner_both_directions_tz_sorted.Symmetrize(second_direction);
                runner_both_directions_tz_sorted.TzSort();
                expect_final_vectors_equivalence(runner_simple, runner_both_directions_tz_sorted);
            }
            // TZ_SORTER + SYMMETRIZER + SYMMETRIZER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                runner::Runner runner_both_directions_tz_sorted_eliminated(model);
                runner_both_directions_tz_sorted_eliminated.TzSort();
                runner_both_directions_tz_sorted_eliminated.Symmetrize(first_direction);
                runner_both_directions_tz_sorted_eliminated.Symmetrize(second_direction);
                runner_both_directions_tz_sorted_eliminated.EliminatePositiveProjections();
                expect_final_vectors_equivalence(
                    runner_simple,
                    runner_both_directions_tz_sorted_eliminated);
            }
            // SYMMETRIZER + SYMMETRIZER + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                runner::Runner runner_both_directions_tz_sorted_eliminated(model);
                runner_both_directions_tz_sorted_eliminated.Symmetrize(first_direction);
                runner_both_directions_tz_sorted_eliminated.Symmetrize(second_direction);
                runner_both_directions_tz_sorted_eliminated.TzSort();
                runner_both_directions_tz_sorted_eliminated.EliminatePositiveProjections();
                expect_final_vectors_equivalence(
                    runner_simple,
                    runner_both_directions_tz_sorted_eliminated);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_three_centers_exchange_triangle(model::Model& model, double first) {
    auto Jfirst = model.getSymbols().addSymbol("J1", first);
    model.getSymbols()
        .assignSymbolToIsotropicExchange(Jfirst, 0, 1)
        .assignSymbolToIsotropicExchange(Jfirst, 1, 2)
        .assignSymbolToIsotropicExchange(Jfirst, 2, 0);

    model.InitializeSSquared();
}

TEST(spectrum_final_equivalence, triangle) {
    std::vector<std::vector<int>> multss = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};
    std::vector<double> js = {10, 17.17, 33};

    for (const auto& mults : multss) {
        for (auto Jfirst : js) {
            group::Group triangle(group::Group::S3, {{2, 0, 1}, {0, 2, 1}});
            model::Model model(mults);
            initialize_three_centers_exchange_triangle(model, Jfirst);

            runner::Runner runner_simple(model);
            runner_simple.BuildSpectra();

            // TZ_SORTER
            {
                runner::Runner runner_tz_sorted(model);
                runner_tz_sorted.TzSort();
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                runner::Runner runner_tz_sorted_eliminated(model);
                runner_tz_sorted_eliminated.TzSort();
                runner_tz_sorted_eliminated.EliminatePositiveProjections();
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted_eliminated);
            }
            // SYMMETRIZER
            {
                runner::Runner runner_symmetrized(model);
                runner_symmetrized.Symmetrize(triangle);
                expect_final_vectors_equivalence(runner_simple, runner_symmetrized);
            }
            {
                group::Group triangle_diff(group::Group::S3, {{1, 2, 0}, {2, 1, 0}});
                runner::Runner runner_symmetrized(model);
                runner_symmetrized.Symmetrize(triangle_diff);
                expect_final_vectors_equivalence(runner_simple, runner_symmetrized);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                runner::Runner runner_symmetrized_tz_sorted(model);
                runner_symmetrized_tz_sorted.TzSort();
                runner_symmetrized_tz_sorted.Symmetrize(triangle);
                expect_final_vectors_equivalence(runner_simple, runner_symmetrized_tz_sorted);
            }
            // SYMMETRIZER + TZ_SORTER
            {
                runner::Runner runner_symmetrized_tz_sorted(model);
                runner_symmetrized_tz_sorted.Symmetrize(triangle);
                runner_symmetrized_tz_sorted.TzSort();
                expect_final_vectors_equivalence(runner_simple, runner_symmetrized_tz_sorted);
            }
            // TZ_SORTER + SYMMETRIZER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                runner::Runner runner_symmetrized_tz_sorted_eliminated(model);
                runner_symmetrized_tz_sorted_eliminated.TzSort();
                runner_symmetrized_tz_sorted_eliminated.Symmetrize(triangle);
                runner_symmetrized_tz_sorted_eliminated.EliminatePositiveProjections();
                expect_final_vectors_equivalence(
                    runner_simple,
                    runner_symmetrized_tz_sorted_eliminated);
            }
            // SYMMETRIZER + TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                runner::Runner runner_symmetrized_tz_sorted_eliminated(model);
                runner_symmetrized_tz_sorted_eliminated.Symmetrize(triangle);
                runner_symmetrized_tz_sorted_eliminated.TzSort();
                runner_symmetrized_tz_sorted_eliminated.EliminatePositiveProjections();
                expect_final_vectors_equivalence(
                    runner_simple,
                    runner_symmetrized_tz_sorted_eliminated);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                runner::Runner runner_symmetrized_tz_sorted_eliminated(model);
                runner_symmetrized_tz_sorted_eliminated.TzSort();
                runner_symmetrized_tz_sorted_eliminated.EliminatePositiveProjections();
                runner_symmetrized_tz_sorted_eliminated.Symmetrize(triangle);
                expect_final_vectors_equivalence(
                    runner_simple,
                    runner_symmetrized_tz_sorted_eliminated);
            }
            //            // SYMMETRIZER + NON_ABELIAN_SIMPLIFIER
            //            {
            //                runner::Runner runner_symmetrized_simplified(mults);
            //                runner_symmetrized_simplified.Symmetrize(triangle);
            //                runner_symmetrized_simplified.NonAbelianSimplify();
            //                initialize_three_centers_exchange_triangle(runner_symmetrized_simplified, Jfirst);
            //                expect_final_vectors_equivalence(runner_simple, runner_symmetrized_simplified);
            //            }
            //            // TZ_SORTER + SYMMETRIZER + NON_ABELIAN_SIMPLIFIER
            //            {
            //                runner::Runner runner_tz_sorted_symmetrized_simplified(mults);
            //                runner_tz_sorted_symmetrized_simplified.TzSort();
            //                runner_tz_sorted_symmetrized_simplified.Symmetrize(triangle);
            //                runner_tz_sorted_symmetrized_simplified.NonAbelianSimplify();
            //                initialize_three_centers_exchange_triangle(runner_tz_sorted_symmetrized_simplified, Jfirst);
            //                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted_symmetrized_simplified);
            //            }
            //            // SYMMETRIZER + TZ_SORTER + NON_ABELIAN_SIMPLIFIER
            //            {
            //                runner::Runner runner_tz_sorted_symmetrized_simplified(mults);
            //                runner_tz_sorted_symmetrized_simplified.Symmetrize(triangle);
            //                runner_tz_sorted_symmetrized_simplified.TzSort();
            //                runner_tz_sorted_symmetrized_simplified.NonAbelianSimplify();
            //                initialize_three_centers_exchange_triangle(runner_tz_sorted_symmetrized_simplified, Jfirst);
            //                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted_symmetrized_simplified);
            //            }
            //            // SYMMETRIZER + NON_ABELIAN_SIMPLIFIER + TZ_SORTER
            //            {
            //                runner::Runner runner_tz_sorted_symmetrized_simplified(mults);
            //                runner_tz_sorted_symmetrized_simplified.Symmetrize(triangle);
            //                runner_tz_sorted_symmetrized_simplified.NonAbelianSimplify();
            //                runner_tz_sorted_symmetrized_simplified.TzSort();
            //                initialize_three_centers_exchange_triangle(runner_tz_sorted_symmetrized_simplified, Jfirst);
            //                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted_symmetrized_simplified);
            //            }
        }
    }
}