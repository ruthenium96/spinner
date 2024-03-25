#include "Parser.h"
#include "Tools.h"

#include <yaml-cpp/yaml.h>

namespace input {

Parser::Parser(const std::string& filename) {
    auto input_node = YAML::LoadFile(filename);
    model_input_parser_.emplace(extractValue<YAML::Node>(input_node, "model_input"));

    optimizations_list_parser_.emplace(extractValue<YAML::Node>(input_node, "optimizations"));

    job_parser_.emplace(extractValue<YAML::Node>(input_node, "job"));

    throw_if_node_is_not_empty(input_node);
}

const std::vector<model::ModelInput>& Parser::getModelInputs() const {
    return model_input_parser_.value().getModelInputs();
}

const std::optional<common::physical_optimization::OptimizationList>&
Parser::getOptimizationList() const {
    return optimizations_list_parser_->getOptimizationList();
}

const std::optional<std::vector<double>>& Parser::getTemperaturesForSimulation() const {
    return job_parser_->getTemperaturesForSimulation();
}

const std::optional<std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>>&
Parser::getNonlinearSolver() const {
    return job_parser_->getNonlinearSolver();
}

const std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>&
Parser::getExperimentalValuesWorker() const {
    return job_parser_->getExperimentalValuesWorker();
}
}  // namespace input