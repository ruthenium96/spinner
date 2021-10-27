#ifndef JULY_RUNNER_H
#define JULY_RUNNER_H

#include <utility>

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
    void BuildMatrix();

    // SPECTRUM OPERATIONS
    void BuildSpectrum();

    [[nodiscard]] const Space& getSpace() const;
    [[nodiscard]] const Matrix& getHamiltonianMatrix() const;
    [[nodiscard]] const std::vector<Matrix>& getNonHamiltonianMatrices() const;
    [[nodiscard]] const Spectrum& getSpectrum() const;

    [[nodiscard]] uint32_t getTotalSpaceSize() const;

  private:
    struct HamiltonianHistory {
        bool has_isotropic_exchange_interactions = false;
    };
    struct SpaceHistory {
        std::vector<Group> applied_groups;
        uint32_t number_of_non_simplified_abelian_groups = 0;
        bool isTzSorted = false;
        bool isNormalized = false; // actually true, if we do not use Symmetrizer
    };

    const spaces::LexicographicIndexConverter converter_;

    Space space_;
    SpaceHistory space_history_;

    Operator hamiltonian_operator_;
    HamiltonianHistory hamiltonian_history_;
    Matrix hamiltonian_matrix_;

    std::vector<Operator> non_hamiltonian_operators;
    std::vector<Matrix> non_hamiltonian_matrices;

    Spectrum spectrum_;
};
} // namespace runner

#endif //JULY_RUNNER_H
