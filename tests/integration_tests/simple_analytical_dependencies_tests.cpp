#include <random>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

double sum_of_s_squared(const lexicographic::IndexConverter& indexConverter) {
    double sum_of_s_squared = 0;
    for (auto spin : indexConverter.get_spins()) {
        sum_of_s_squared += spin * (spin + 1);
    }
    return sum_of_s_squared;
}

TEST(simple_analytical_dependencies, nothing) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> Theta_dist(-100, -1.0);

    std::vector<std::vector<int>> mults_cases =
        {{2}, {3}, {4}, {5}, {2, 2}, {2, 2, 2, 2}, {2, 3, 4, 5}};

    for (const auto& mults : mults_cases) {
        for (size_t _ = 0; _ < 10; ++_) {
            const double g_exact = 2.0;

            std::vector<magnetic_susceptibility::ValueAtTemperature> values;

            {
                model::Model model(mults);
                double g_value = g_exact;
                auto g = model.getSymbols().addSymbol("g", g_value);
                for (size_t i = 0; i < mults.size(); ++i) {
                    model.getSymbols().assignSymbolToGFactor(g, i);
                }
                model.InitializeSSquared();

                runner::Runner runner(model);

                runner.BuildSpectra();
                runner.BuildMuSquaredWorker();

                double sum_of_s_squared_ = sum_of_s_squared(runner.getIndexConverter());

                for (size_t temperature = 1; temperature < 301; ++temperature) {
                    double calculated_value =
                        runner.getMuSquaredWorker().calculateTheoreticalMuSquared(temperature);
                    double exact_value = (g_exact * g_exact * sum_of_s_squared_);
                    // TODO: epsilon
                    EXPECT_NEAR(calculated_value, exact_value, 1e-9);
                }
            }
        }
    }
}

TEST(simple_analytical_dependencies, Theta) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> Theta_dist(-100, -1.0);

    std::vector<std::vector<int>> mults_cases =
        {{2}, {3}, {4}, {5}, {2, 2}, {2, 2, 2, 2}, {2, 3, 4, 5}};

    for (const auto& mults : mults_cases) {
        for (size_t _ = 0; _ < 10; ++_) {
            const double Theta_exact = Theta_dist(rng);
            const double g_exact = 2.0;

            std::vector<magnetic_susceptibility::ValueAtTemperature> values;

            {
                model::Model model(mults);
                auto Theta = model.getSymbols().addSymbol("Theta", Theta_exact);
                model.getSymbols().assignSymbolToTheta(Theta);
                double g_value = g_exact;
                auto g = model.getSymbols().addSymbol("g", g_value);
                for (size_t i = 0; i < mults.size(); ++i) {
                    model.getSymbols().assignSymbolToGFactor(g, i);
                }
                model.InitializeSSquared();

                runner::Runner runner(model);

                runner.BuildSpectra();
                runner.BuildMuSquaredWorker();

                double sum_of_s_squared_ = sum_of_s_squared(runner.getIndexConverter());

                for (size_t temperature = 1; temperature < 301; ++temperature) {
                    double calculated_value =
                        runner.getMuSquaredWorker().calculateTheoreticalMuSquared(temperature);
                    double exact_value = (g_exact * g_exact * sum_of_s_squared_) * temperature
                        / (temperature - Theta_exact);
                    // TODO: epsilon
                    EXPECT_NEAR(calculated_value, exact_value, 1e-9);
                }
            }
        }
    }
}