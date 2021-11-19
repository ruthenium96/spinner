#ifndef JULY_RUNNER_H
#define JULY_RUNNER_H

#include <utility>

#include "common/Quantity.h"
#include "common/symbols/Symbols.h"
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
    void AddIsotropicExchange(double initial_value, size_t center_a, size_t center_b);

    // OPERATOR OPERATIONS
    void InitializeSSquared();
    void FinalizeIsotropicInteraction();

    // MATRIX OPERATIONS
    void BuildMatrices();

    // SPECTRUM OPERATIONS
    void BuildSpectra();

    [[nodiscard]] const lexicographic::IndexConverter& getIndexConverter() const;
    [[nodiscard]] const Operator& getOperator(common::QuantityEnum) const;
    [[nodiscard]] const Space& getSpace() const;
    [[nodiscard]] const Spectrum& getSpectrum(common::QuantityEnum) const;
    [[nodiscard]] const Matrix& getMatrix(common::QuantityEnum) const;

    [[nodiscard]] uint32_t getTotalSpaceSize() const;

  private:
    struct MatrixHistory {
        bool matrices_was_built = false;
    };
    struct HamiltonianOperatorHistory {
        bool has_isotropic_exchange_interactions_initialized = false;
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

    // TODO: should we split these containers?
    std::map<common::QuantityEnum, Operator> operators_;
    std::map<common::QuantityEnum, Matrix> matrices_;
    std::map<common::QuantityEnum, Spectrum> spectra_;

    void BuildSpectraUsingMatrices();
    void BuildSpectraWithoutMatrices();

    MatrixHistory matrix_history_;
    HamiltonianOperatorHistory hamiltonian_history_;
    SpaceHistory space_history_;
};
}  // namespace runner

#endif  //JULY_RUNNER_H
