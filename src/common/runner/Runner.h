#ifndef JULY_RUNNER_H
#define JULY_RUNNER_H

#include <utility>

#include "ConsistentModelOptimizationList.h"
#include "src/common/Quantity.h"
#include "src/entities/magnetic_susceptibility/MuSquaredWorker.h"
#include "src/entities/magnetic_susceptibility/UniqueGOnlySSquaredMuSquaredWorker.h"
#include "src/entities/matrix/Matrix.h"
#include "src/entities/spectrum/Spectrum.h"
#include "src/space/Space.h"

namespace runner {
class Runner {
  public:
    Runner(model::Model model, common::physical_optimization::OptimizationList optimizationList);
    // constructor for Runner with no optimizations:
    explicit Runner(model::Model model);

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
    void minimizeResidualError();

    const lexicographic::IndexConverter& getIndexConverter() const;
    const model::operators::Operator& getOperator(common::QuantityEnum) const;
    const space::Space& getSpace() const;
    const Spectrum& getSpectrum(common::QuantityEnum) const;
    const Matrix& getMatrix(common::QuantityEnum) const;
    const model::operators::Operator& getOperatorDerivative(
        common::QuantityEnum,
        model::symbols::SymbolTypeEnum,
        const model::symbols::SymbolName&) const;
    const Spectrum& getSpectrumDerivative(
        common::QuantityEnum,
        model::symbols::SymbolTypeEnum,
        const model::symbols::SymbolName&) const;
    const Matrix& getMatrixDerivative(
        common::QuantityEnum,
        model::symbols::SymbolTypeEnum,
        const model::symbols::SymbolName&) const;
    const magnetic_susceptibility::MuSquaredWorker& getMuSquaredWorker() const;
    const model::symbols::Symbols& getSymbols() const;

  private:
    ConsistentModelOptimizationList consistentModelOptimizationList_;
    const space::Space space_;

    const model::Model& getModel() const;
    model::Model& getModel();
    const common::physical_optimization::OptimizationList& getOptimizationList() const;

    struct MatrixHistory {
        bool matrices_was_built = false;
    };

    common::Quantity energy;
    std::optional<common::Quantity> s_squared;
    std::map<model::symbols::SymbolName, common::Quantity>
        derivative_of_energy_wrt_exchange_parameters;

    void stepOfRegression(
        const std::vector<model::symbols::SymbolName>&,
        const std::vector<double>&,
        double&,
        std::vector<double>&);

    // TODO: can we use std::optional<magnetic_susceptibility::MuSquaredWorker> instead?
    std::optional<std::unique_ptr<magnetic_susceptibility::MuSquaredWorker>> mu_squared_worker;
    std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>
        experimental_values_worker_;

    void BuildSpectraUsingMatrices(size_t number_of_blocks);
    void BuildSpectraWithoutMatrices(size_t number_of_blocks);

    MatrixHistory matrix_history_;
};
}  // namespace runner

#endif  //JULY_RUNNER_H
