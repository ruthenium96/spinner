#include "Executer.h"

namespace runner {
void Executer::execute(input::Parser parser) {
    auto optimization_list = parser.getOptimizationList().value();

    for (const auto& model_input : parser.getModelInputs()) {
        auto runner = runner::Runner(model_input, optimization_list);

        if (parser.getTemperaturesForSimulation().has_value()) {
            for (auto temp : parser.getTemperaturesForSimulation().value()) {
                std::cout << temp << " " <<
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(temp)
                          << std::endl;
            }
        } else if (parser.getNonlinearSolver().has_value()) {
            runner.initializeExperimentalValues(parser.getExperimentalValuesWorker().value());

            runner.minimizeResidualError(parser.getNonlinearSolver().value());

            auto theor_final = runner.getMagneticSusceptibilityController().getTheoreticalValues();

            std::cout << std::endl;
            for (auto [temperature, value] : theor_final) {
                std::cout << temperature << " " << value << std::endl;
            }
        }
    }
}
}  // namespace runner