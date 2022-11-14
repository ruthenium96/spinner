#include "abstractSolver_tests.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::optimNMAdapter>() {
    return std::make_shared<nonlinear_solver::optimNMAdapter>();
};

typedef testing::Types<nonlinear_solver::optimNMAdapter> optimNM;
INSTANTIATE_TYPED_TEST_SUITE_P(optimNMSolverTests, fitting_magnetic_susceptibility_simple, optimNM);

TEST(fitting_magnetic_susceptibility_advanced, optimNM_fit_theoretical_curve_2222_JAF_g_g) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-50, -40);
    std::uniform_real_distribution<double> g_one_dist(1.8, 2.0);
    std::uniform_real_distribution<double> g_two_dist(2.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double g_one_exact = g_one_dist(rng);
        const double g_two_exact = g_two_dist(rng);
        std::vector<int> mults = {2, 2, 2, 2};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            model::ModelInput model(mults);
            double J_value = J_exact;
            auto J = model.modifySymbolicWorker().addSymbol("J", J_value);
            model.modifySymbolicWorker()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);
            double g_one_value = g_one_exact;
            double g_two_value = g_two_exact;
            auto g_one = model.modifySymbolicWorker().addSymbol("g1", g_one_value);
            auto g_two = model.modifySymbolicWorker().addSymbol("g2", g_two_value);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.modifySymbolicWorker().assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.modifySymbolicWorker().assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            model::ModelInput model(mults);
            double J_value = -45.0;
            auto J = model.modifySymbolicWorker().addSymbol("J", J_value);
            model.modifySymbolicWorker()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);
            double g_one_value = 1.8;
            double g_two_value = 3.0;
            auto g_one = model.modifySymbolicWorker().addSymbol("g1", g_one_value);
            auto g_two = model.modifySymbolicWorker().addSymbol("g2", g_two_value);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.modifySymbolicWorker().assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.modifySymbolicWorker().assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(std::make_shared<nonlinear_solver::optimNMAdapter>());

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double J_fitted = runner.getSymbolicWorker().getValueOfName(J);
            double g_one_fitted = runner.getSymbolicWorker().getValueOfName(g_one);
            double g_two_fitted = runner.getSymbolicWorker().getValueOfName(g_two);
            double J_range = std::abs(J_fitted / 1000);
            double g_one_range = std::abs((g_one_fitted) / 1000);
            double g_two_range = std::abs((g_two_fitted) / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(J_fitted, J_exact, J_range);
            EXPECT_NEAR(g_one_fitted, g_one_exact, g_one_range);
            EXPECT_NEAR(g_two_fitted, g_two_exact, g_two_range);
        }
    }
}

TEST(
    fitting_magnetic_susceptibility_advanced,
    optimNM_fit_theoretical_curve_2222_JAF_g_g_unique_vs_different) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-50, -40);
    std::uniform_real_distribution<double> g_dist(1.8, 2.2);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double g_exact = g_dist(rng);
        std::vector<int> mults = {2, 2, 2, 2};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            model::ModelInput model(mults);
            double J_value = J_exact;
            auto J = model.modifySymbolicWorker().addSymbol("J", J_value);
            model.modifySymbolicWorker()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);
            double g_value = g_exact;
            auto g = model.modifySymbolicWorker().addSymbol("g1", g_value);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.modifySymbolicWorker().assignSymbolToGFactor(g, i);
            }

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            model::ModelInput model(mults);
            double J_value = -45.0;
            auto J = model.modifySymbolicWorker().addSymbol("J", J_value);
            model.modifySymbolicWorker()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);
            double g_one_value = 1.8;
            double g_two_value = 2.2;
            auto g_one = model.modifySymbolicWorker().addSymbol("g1", g_one_value);
            auto g_two = model.modifySymbolicWorker().addSymbol("g2", g_two_value);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.modifySymbolicWorker().assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.modifySymbolicWorker().assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(std::make_shared<nonlinear_solver::optimNMAdapter>());

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double J_fitted = runner.getSymbolicWorker().getValueOfName(J);
            double g_one_fitted = runner.getSymbolicWorker().getValueOfName(g_one);
            double g_two_fitted = runner.getSymbolicWorker().getValueOfName(g_two);
            double J_range = std::abs(J_fitted / 1000);
            double g_one_range = std::abs((g_one_fitted) / 1000);
            double g_two_range = std::abs((g_two_fitted) / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(J_fitted, J_exact, J_range);
            EXPECT_NEAR(g_one_fitted, g_exact, g_one_range);
            EXPECT_NEAR(g_two_fitted, g_exact, g_two_range);
        }
    }
}
