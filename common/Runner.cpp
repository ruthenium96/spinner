#include "Runner.h"
#include <components/matrix/MatrixBuilder.h>
#include "components/operator/IsotropicExchangeHamiltonian.h"
#include "components/space/NonAbelianSimplifier.h"
#include "components/space/Symmetrizer.h"
#include "components/space/TzSorter.h"


runner::Runner::Runner(std::vector<int> mults) : converter_(std::move(mults)), space_(converter_.total_space_size)
{}

void runner::Runner::NonAbelianSimplify() {
    if (space_history_.number_of_non_simplified_abelian_groups == 0) {
        return;
    }
    if (space_history_.number_of_non_simplified_abelian_groups != 1) {
        throw std::invalid_argument("Non-Abelian simplification after using of two Non-Abelian Symmetrizers "
                                    "currently is not allowed. Use Non-Abelian simplification twice.");
    }
    NonAbelianSimplifier nonAbelianSimplifier;
    space_ = nonAbelianSimplifier.apply(std::move(space_));
    space_history_.number_of_non_simplified_abelian_groups = 0;
}

void runner::Runner::Symmetrize(Group new_group) {
    // check if user trying to use the same Group for a second time:
    if (std::count(space_history_.applied_groups.begin(), space_history_.applied_groups.end(), new_group)) {
        return;
    }

    Symmetrizer symmetrizer(converter_, new_group);
    space_ = symmetrizer.apply(std::move(space_));

    if (!new_group.properties.is_abelian) {
        ++space_history_.number_of_non_simplified_abelian_groups;
    }
    space_history_.applied_groups.emplace_back(std::move(new_group));
}

void runner::Runner::Symmetrize(Group::GroupTypeEnum group_name, std::vector<Permutation> generators) {
    Group new_group(group_name, std::move(generators));
    Symmetrize(new_group);
}

void runner::Runner::TzSort() {
    // It does not make any sense to use tz_sorter twice.
    if (space_history_.isTzSorted) {
        return;
    }
    TzSorter tz_sorter(converter_);
    space_ = tz_sorter.apply(std::move(space_));
    space_history_.isTzSorted = true;
}

const Space &runner::Runner::getSpace() const {
    return space_;
}

uint32_t runner::Runner::getTotalSpaceSize() const {
    return converter_.total_space_size;
}

void runner::Runner::AddIsotropicExchange(arma::dmat isotropic_exchange_parameters) {
    if (hamiltonian_history_.has_isotropic_exchange_interactions) {
        throw std::invalid_argument("Trying to add isotropic exchange twice");
    }

    hamiltonian_.two_center_terms.emplace_back(new IsotropicExchangeHamiltonian(std::move(isotropic_exchange_parameters)));

    hamiltonian_history_.has_isotropic_exchange_interactions = true;
}

void runner::Runner::BuildMatrix() {
    if (!space_history_.isNormalized) {
        for (auto& subspace : space_.blocks) {
            // TODO: maybe, we can implement normalize as Space method
            subspace.decomposition.normalize();
        }
        space_history_.isNormalized = true;
    }

    MatrixBuilder matrix_builder(converter_);
    matrix_builder.apply(space_, hamiltonian_);
}
