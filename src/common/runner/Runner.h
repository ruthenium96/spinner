#ifndef JULY_RUNNER_H
#define JULY_RUNNER_H

#include <utility>

#include "src/common/Quantity.h"
#include "src/entities/magnetic_susceptibility/MuSquaredWorker.h"
#include "src/entities/magnetic_susceptibility/UniqueGOnlySSquaredMuSquaredWorker.h"
#include "src/entities/matrix/Matrix.h"
#include "src/entities/space/Space.h"
#include "src/entities/spectrum/Spectrum.h"
#include "src/group/Group.h"
#include "src/model/Model.h"

// TODO: check symmetry consistency of symbols and groups.
namespace runner {
class Runner {
  public:
    explicit Runner(model::Model model);

    //    explicit Runner(const std::vector<int>& mults);

    // SPACE OPERATIONS
    void EliminatePositiveProjections();

    void NonAbelianSimplify();

    void Symmetrize(group::Group new_group);
    void
    Symmetrize(group::Group::GroupTypeEnum group_name, std::vector<group::Permutation> generators);

    void TzSort();

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
    const Space& getSpace() const;
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
    bool model_is_finished = false;

    model::Model model_;

    struct MatrixHistory {
        bool matrices_was_built = false;
    };
    struct SpaceHistory {
        std::vector<group::Group> applied_groups;
        uint32_t number_of_non_simplified_abelian_groups = 0;
        bool isTzSorted = false;
        bool isNormalized = false;  // actually true, if we do not use Symmetrizer
        bool isNonAbelianSimplified = false;
        bool isPositiveProjectionsEliminated = false;
    };

    Space space_;

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

    // TODO: delete these functions after refactoring.
    void finish_the_model();
    void throw_if_model_is_finished(const std::string& error);

    MatrixHistory matrix_history_;
    SpaceHistory space_history_;
};
}  // namespace runner

#endif  //JULY_RUNNER_H
