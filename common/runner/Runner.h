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
    void NonAbelianSimplify();

    void Symmetrize(group::Group new_group);
    void
    Symmetrize(group::Group::GroupTypeEnum group_name, std::vector<group::Permutation> generators);

    void TzSort();

    // SYMBOLS OPERATIONS
    symbols::SymbolName AddSymbol(
        const std::string& name,
        double initial_value,
        bool is_changeable,
        symbols::SymbolTypeEnum type_enum);
    symbols::SymbolName
    AddSymbol(const std::string& name, double initial_value, bool is_changeable);
    symbols::SymbolName AddSymbol(const std::string& name, double initial_value);
    // TODO: Symbol& modifySymbol();

    // OPERATOR OPERATIONS
    void AssignSymbolToIsotropicExchange(
        const symbols::SymbolName& symbol_name,
        size_t center_a,
        size_t center_b);
    void AssignSymbolToGFactor(
        const symbols::SymbolName& symbol_name,
        size_t center_a);  // actually, it is not "Operator" operation
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

    [[nodiscard]] const lexicographic::IndexConverter& getIndexConverter() const;
    [[nodiscard]] const Operator& getOperator(common::QuantityEnum) const;
    [[nodiscard]] const Space& getSpace() const;
    [[nodiscard]] const Spectrum& getSpectrum(common::QuantityEnum) const;
    [[nodiscard]] const Matrix& getMatrix(common::QuantityEnum) const;
    [[nodiscard]] const Operator& getOperatorDerivative(
        common::QuantityEnum,
        symbols::SymbolTypeEnum,
        const symbols::SymbolName&) const;
    [[nodiscard]] const Spectrum& getSpectrumDerivative(
        common::QuantityEnum,
        symbols::SymbolTypeEnum,
        const symbols::SymbolName&) const;
    [[nodiscard]] const Matrix& getMatrixDerivative(
        common::QuantityEnum,
        symbols::SymbolTypeEnum,
        const symbols::SymbolName&) const;
    [[nodiscard]] const std::unique_ptr<magnetic_susceptibility::MuSquaredWorker>&
    getPtrToMuSquaredWorker() const;
    [[nodiscard]] const symbols::Symbols& getSymbols() const;

  private:
    struct MatrixHistory {
        bool matrices_was_built = false;
    };
    struct HamiltonianOperatorHistory {
        bool has_isotropic_exchange_interactions_initialized = false;
        bool has_isotropic_exchange_interactions_finalized = false;
    };
    struct SpaceHistory {
        std::vector<group::Group> applied_groups;
        uint32_t number_of_non_simplified_abelian_groups = 0;
        bool isTzSorted = false;
        bool isNormalized = false;  // actually true, if we do not use Symmetrizer
        bool isNonAbelianSimplified = false;
    };

    const lexicographic::IndexConverter converter_;
    symbols::Symbols symbols_;

    Space space_;

    Operator operator_energy;
    Operator operator_s_squared;
    std::map<symbols::SymbolName, Operator> operator_derivative_of_energy_wrt_exchange_parameters;
    Matrix matrix_energy;
    Matrix matrix_s_squared;
    std::map<symbols::SymbolName, Matrix> matrix_derivative_of_energy_wrt_exchange_parameters;
    Spectrum spectrum_energy;
    Spectrum spectrum_s_squared;
    std::map<symbols::SymbolName, Spectrum> spectrum_derivative_of_energy_wrt_exchange_parameters;

    // TODO: can we use std::optional<magnetic_susceptibility::MuSquaredWorker> instead?
    std::optional<std::unique_ptr<magnetic_susceptibility::MuSquaredWorker>> mu_squared_worker;
    std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>
        experimental_values_worker_;

    void BuildSpectraUsingMatrices();
    void BuildSpectraWithoutMatrices();

    MatrixHistory matrix_history_;
    HamiltonianOperatorHistory hamiltonian_history_;
    SpaceHistory space_history_;
};
}  // namespace runner

#endif  //JULY_RUNNER_H
