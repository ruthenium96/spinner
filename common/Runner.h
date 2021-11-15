#ifndef JULY_RUNNER_H
#define JULY_RUNNER_H

#include <utility>

#include "common/Quantity.h"
#include "entities/matrix/Matrix.h"
#include "entities/operator/Operator.h"
#include "entities/space/Space.h"
#include "entities/spectrum/Spectrum.h"
#include "groups/Group.h"

namespace runner {
class Runner {
  public:
    explicit Runner(std::vector<int> mults);

    // SPACE OPERATIONS
    void NonAbelianSimplify();

    void Symmetrize(Group new_group);
    void Symmetrize(Group::GroupTypeEnum group_name, std::vector<Permutation> generators);

    void TzSort();

    // OPERATOR OPERATIONS
    void AddIsotropicExchange(arma::dmat isotropic_exchange_parameters);
    void InitializeSSquared();

    // MATRIX OPERATIONS
    void BuildMatrices();

    // SPECTRUM OPERATIONS
    void BuildSpectra();

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
        bool has_isotropic_exchange_interactions = false;
    };
    struct SpaceHistory {
        std::vector<Group> applied_groups;
        uint32_t number_of_non_simplified_abelian_groups = 0;
        bool isTzSorted = false;
        bool isNormalized = false; // actually true, if we do not use Symmetrizer
        bool isNonAbelianSimplified = false;
    };

    const spaces::LexicographicIndexConverter converter_;

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
} // namespace runner

#endif //JULY_RUNNER_H
