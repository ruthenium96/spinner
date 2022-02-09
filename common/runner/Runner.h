#ifndef JULY_RUNNER_H
#define JULY_RUNNER_H

#include <utility>

#include "common/Quantity.h"
#include "common/symbols/Symbols.h"
#include "entities/magnetic_susceptibility/MuSquaredWorker.h"
#include "entities/magnetic_susceptibility/UniqueGOnlySSquaredMuSquaredWorker.h"
#include "entities/matrix/Matrix.h"
#include "entities/operator/Operator.h"
#include "entities/space/Space.h"
#include "entities/spectrum/Spectrum.h"
#include "group/Group.h"

namespace runner {
class Runner {
  public:
    explicit Runner(const std::vector<int>& mults);

    // SPACE OPERATIONS
    void EliminatePositiveProjections();

    void NonAbelianSimplify();

    void Symmetrize(group::Group new_group);
    void
    Symmetrize(group::Group::GroupTypeEnum group_name, std::vector<group::Permutation> generators);

    void TzSort();

    // SYMBOLS OPERATIONS

    // OPERATOR OPERATIONS
    void InitializeIsotropicExchange();
    void InitializeIsotropicExchangeDerivatives();
    void InitializeSSquared();

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
    std::map<symbols::SymbolName, double> calculateTotalDerivatives();
    void minimizeResidualError();

    const lexicographic::IndexConverter& getIndexConverter() const;
    const Operator& getOperator(common::QuantityEnum) const;
    const Space& getSpace() const;
    const Spectrum& getSpectrum(common::QuantityEnum) const;
    const Matrix& getMatrix(common::QuantityEnum) const;
    const Operator& getOperatorDerivative(
        common::QuantityEnum,
        symbols::SymbolTypeEnum,
        const symbols::SymbolName&) const;
    const Spectrum& getSpectrumDerivative(
        common::QuantityEnum,
        symbols::SymbolTypeEnum,
        const symbols::SymbolName&) const;
    const Matrix& getMatrixDerivative(
        common::QuantityEnum,
        symbols::SymbolTypeEnum,
        const symbols::SymbolName&) const;
    const magnetic_susceptibility::MuSquaredWorker& getMuSquaredWorker() const;
    const symbols::Symbols& getConstSymbols() const;

    symbols::Symbols& getMutableSymbols();

  private:
    bool model_is_finished = false;
    struct MatrixHistory {
        bool matrices_was_built = false;
    };
    struct OperatorsHistory {
        bool isotropic_exchange_in_hamiltonian = false;
        bool s_squared = false;
        bool isotropic_exchange_derivatives = false;
    };
    struct SpaceHistory {
        std::vector<group::Group> applied_groups;
        uint32_t number_of_non_simplified_abelian_groups = 0;
        bool isTzSorted = false;
        bool isNormalized = false;  // actually true, if we do not use Symmetrizer
        bool isNonAbelianSimplified = false;
        bool isPositiveProjectionsEliminated = false;
    };

    const lexicographic::IndexConverter converter_;
    symbols::Symbols symbols_;

    Space space_;

    common::Quantity energy;
    std::optional<common::Quantity> s_squared;
    std::map<symbols::SymbolName, common::Quantity> derivative_of_energy_wrt_exchange_parameters;
    
    void stepOfRegression(
        const std::vector<symbols::SymbolName>&,
        const std::vector<double>&,
        double&,
        std::vector<double>&);

    // TODO: can we use std::optional<magnetic_susceptibility::MuSquaredWorker> instead?
    std::optional<std::unique_ptr<magnetic_susceptibility::MuSquaredWorker>> mu_squared_worker;
    std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>
        experimental_values_worker_;

    void BuildSpectraUsingMatrices(size_t number_of_blocks);
    void BuildSpectraWithoutMatrices(size_t number_of_blocks);

    void finish_the_model();
    void throw_if_model_is_finished(const std::string& error);

    MatrixHistory matrix_history_;
    OperatorsHistory operators_history_;
    SpaceHistory space_history_;
};
}  // namespace runner

#endif  //JULY_RUNNER_H
