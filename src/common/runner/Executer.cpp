#include "Executer.h"

#include "src/common/Logger.h"

namespace runner {
void Executer::execute(input::Parser parser) {
    auto optimization_list = parser.getOptimizationList().value();

    for (const auto& model_input : parser.getModelInputs()) {
        auto runner = runner::Runner(model_input, optimization_list);

        if (parser.getTemperaturesForSimulation().has_value()) {
            for (auto temp : parser.getTemperaturesForSimulation().value()) {
                double value = runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(temp);
                common::Logger::basic("{:.8}    {:.8}", temp, value);
            }
        } else if (parser.getNonlinearSolver().has_value()) {
            runner.initializeExperimentalValues(parser.getExperimentalValuesWorker().value());

            runner.minimizeResidualError(parser.getNonlinearSolver().value());

            common::Logger::basic_msg("\nMuSquared, Temperature");
            auto theor_final = runner.getMagneticSusceptibilityController().getTheoreticalValues();
            for (auto [temperature, value] : theor_final) {
                common::Logger::basic("{:.8e}    {:.8e}", temperature, value);
            }
        }
    }
    common::Logger::separate(0, common::basic);
}
}  // namespace runner