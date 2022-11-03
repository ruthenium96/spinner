#ifndef SPINNER_RUNNER_H
#define SPINNER_RUNNER_H

#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>

#include "ConsistentModelOptimizationList.h"
#include "src/common/Quantity.h"
#include "src/entities/data_structures/AbstractFactory.h"
#include "src/entities/magnetic_susceptibility/MagneticSusceptibilityController.h"
#include "src/entities/matrix/Matrix.h"
#include "src/entities/spectrum/Spectrum.h"
#include "src/nonlinear_solver/AbstractNonlinearSolver.h"
#include "src/space/Space.h"

namespace runner {
class Runner {
  public:
    Runner(
        model::ModelInput model,
        common::physical_optimization::OptimizationList optimizationList,
        std::shared_ptr<quantum::linear_algebra::AbstractFactory> algebraDataFactory);
    // constructor for Runner with default algebra package:
    Runner(
        model::ModelInput model,
        common::physical_optimization::OptimizationList optimizationList);
    // constructor for Runner with no optimizations:
    Runner(
        model::ModelInput model,
        std::shared_ptr<quantum::linear_algebra::AbstractFactory> algebraDataFactory);
    // constructor for Runner with no optimizations and default algebra package:
    explicit Runner(model::ModelInput model);

    // TODO: This function is public only for tests. Fix it?
    void initializeDerivatives();

    // SPACE OPERATIONS

    // SYMBOLS OPERATIONS

    // MATRIX OPERATIONS
    void BuildMatrices();

    // SPECTRUM OPERATIONS
    void BuildSpectra();

    // CHIT OPERATIONS
    void BuildMuSquaredWorker();
    void initializeExperimentalValues(
        const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_data,
        magnetic_susceptibility::ExperimentalValuesEnum experimental_quantity_type,
        double number_of_centers_ratio);
    std::map<model::symbols::SymbolName, double> calculateTotalDerivatives();
    void minimizeResidualError(std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>);

    const lexicographic::IndexConverter& getIndexConverter() const;
    const model::operators::Operator& getOperator(common::QuantityEnum) const;
    const space::Space& getSpace() const;
    const Spectrum& getSpectrum(common::QuantityEnum) const;
    const Matrix& getMatrix(common::QuantityEnum) const;
    const model::operators::Operator& getOperatorDerivative(
        common::QuantityEnum,
        const model::symbols::SymbolName&) const;
    const Spectrum& getSpectrumDerivative(
        common::QuantityEnum,
        const model::symbols::SymbolName&) const;
    const Matrix& getMatrixDerivative(
        common::QuantityEnum,
        const model::symbols::SymbolName&) const;
    const magnetic_susceptibility::MagneticSusceptibilityController&
    getMagneticSusceptibilityController() const;
    const model::symbols::Symbols& getSymbols() const;

    std::shared_ptr<quantum::linear_algebra::AbstractFactory> getAlgebraDataFactory() const;

  private:
    ConsistentModelOptimizationList consistentModelOptimizationList_;
    const space::Space space_;
    // TODO: rename it.
    std::shared_ptr<quantum::linear_algebra::AbstractFactory> algebraDataFactory_;

    const model::Model& getModel() const;
    model::Model& getModel();
    const common::physical_optimization::OptimizationList& getOptimizationList() const;

    struct MatrixHistory {
        bool matrices_was_built = false;
    };

    common::Quantity energy;
    std::optional<common::Quantity> s_squared;
    std::optional<common::Quantity> g_sz_squared;
    std::map<std::pair<common::QuantityEnum, model::symbols::SymbolName>, common::Quantity>
        derivatives_map_;

    double stepOfRegression(
        const std::vector<model::symbols::SymbolName>&,
        const std::vector<double>&,
        std::vector<double>&,
        bool isGradientRequired);

    std::optional<magnetic_susceptibility::MagneticSusceptibilityController>
        magnetic_susceptibility_controller_;
    std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>
        experimental_values_worker_;

    void BuildSpectraUsingMatrices(size_t number_of_blocks);
    void BuildSpectraWithoutMatrices(size_t number_of_blocks);

    MatrixHistory matrix_history_;
};
}  // namespace runner

#endif  //SPINNER_RUNNER_H
