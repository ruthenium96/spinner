#include "Runner.h"
#include "components/matrix/MatrixBuilder.h"
#include "components/operator/ConstantOperator.h"
#include "components/operator/ScalarProduct.h"
#include "components/space/NonAbelianSimplifier.h"
#include "components/space/Symmetrizer.h"
#include "components/space/TzSorter.h"
#include "components/spectrum/SpectrumBuilder.h"

namespace {
    // TODO: It requires three allocation instead of two. Rewrite it.
    bool OperatorParametersMatchSymmetries(const std::vector<Group>& applied_groups, const arma::dmat initial_parameters) {
        for (const auto& group : applied_groups) {
            for (const auto& element : group.elements_) {
                arma::dmat parameters_rows_swaped;
                parameters_rows_swaped.resize(arma::size(initial_parameters));
                for (size_t i = 0; i < element.size(); ++i) {
                    parameters_rows_swaped.row(i) = initial_parameters.row(element[i]);
                }
                arma::dmat parameters_all_swaped;
                parameters_all_swaped.resize(arma::size(initial_parameters));
                for (size_t i = 0; i < element.size(); ++i) {
                    parameters_all_swaped.col(i) = parameters_rows_swaped.col(element[i]);
                }
                // TODO: EPSILON
                if (!arma::approx_equal(parameters_all_swaped, initial_parameters, "absdiff", 0.0000001)) {
                    return false;
                }
            }
        }
        return true;
    }
}

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

    // TODO: Should we move it to constructor?
    if (operators_.count(QuantityEnum::Energy) == 0) {
        operators_[QuantityEnum::Energy] = Operator();
    }

    operators_.at(QuantityEnum::Energy).two_center_terms
    .emplace_back(new ScalarProduct(std::move(isotropic_exchange_parameters)));

    hamiltonian_history_.has_isotropic_exchange_interactions = true;
}

void runner::Runner::BuildMatrices() {

    for (const auto& ptr_to_term : operators_.at(QuantityEnum::Energy).two_center_terms) {
        if (!OperatorParametersMatchSymmetries(space_history_.applied_groups,
                                               ptr_to_term->get_parameters().replace(arma::datum::nan, 0))) {
            throw std::invalid_argument("Operator parameters does not match applied symmetries");
        }
    }

    if (!space_history_.isNormalized) {
        for (auto& subspace : space_.blocks) {
            // TODO: maybe, we can implement normalize as Space method
            subspace.decomposition.normalize();
        }
        space_history_.isNormalized = true;
    }

    MatrixBuilder matrix_builder(converter_);
    for (const auto& pair : operators_) {
        matrices_[pair.first] = matrix_builder.apply(space_, pair.second);
    }

    matrix_history_.matrices_was_built = true;
}

void runner::Runner::InitializeSSquared() {
    Operator s_squared_operator_;
    double sum_of_s_squared = 0;
    for (double spin : converter_.get_spins()) {
        sum_of_s_squared += spin * (spin + 1);
    }
    s_squared_operator_.zero_center_terms.emplace_back(new ConstantOperator(sum_of_s_squared));
    s_squared_operator_.two_center_terms.emplace_back(new ScalarProduct(converter_.get_spins().size()));

    operators_[QuantityEnum::S_total_squared] = std::move(s_squared_operator_);
}


void runner::Runner::BuildSpectra() {
    if (matrix_history_.matrices_was_built) {
        BuildSpectraUsingMatrices();
    } else {
        BuildSpectraWithoutMatrices();
    }
}

void runner::Runner::BuildSpectraUsingMatrices() {
    SpectrumBuilder spectrumBuilder;

    size_t number_of_blocks = space_.blocks.size();

    for (const auto& pair : matrices_) {
        spectra_[pair.first] = Spectrum();
        spectra_.at(pair.first).blocks.resize(number_of_blocks);
    }

    for (size_t block = 0; block < number_of_blocks; ++block) {
        DenseMatrix unitary_transformation_matrix;
        spectra_.at(QuantityEnum::Energy).blocks[block] =
                spectrumBuilder.apply_to_subentity_energy(
                        matrices_.at(QuantityEnum::Energy).blocks[block],
                        unitary_transformation_matrix
                );

        for (const auto& pair : matrices_) {
            if (pair.first != QuantityEnum::Energy) {
                spectra_.at(pair.first).blocks[block] =
                        spectrumBuilder.apply_to_subentity_non_energy(
                                pair.second.blocks[block],
                                unitary_transformation_matrix
                                );
            }
        }
    }
}


void runner::Runner::BuildSpectraWithoutMatrices() {
    MatrixBuilder matrixBuilder(converter_);
    SpectrumBuilder spectrumBuilder;

    size_t number_of_blocks = space_.blocks.size();

    for (const auto& pair : operators_) {
        spectra_[pair.first] = Spectrum();
        spectra_.at(pair.first).blocks.resize(number_of_blocks);
    }

    if (!space_history_.isNormalized) {
        for (auto& subspace : space_.blocks) {
            // TODO: maybe, we can implement normalize as Space method
            subspace.decomposition.normalize();
        }
        space_history_.isNormalized = true;
    }

    for (size_t block = 0; block < number_of_blocks; ++block) {
        DenseMatrix unitary_transformation_matrix;

        {
            Submatrix hamiltonian_submatrix = matrixBuilder
                    .apply_to_subentity(space_.blocks[block],
                                        operators_.at(QuantityEnum::Energy));
            spectra_.at(QuantityEnum::Energy).blocks[block] =
                    spectrumBuilder.apply_to_subentity_energy(hamiltonian_submatrix,
                                                              unitary_transformation_matrix);
        }

        for (const auto& pair : operators_) {
            if (pair.first != QuantityEnum::Energy) {
                Submatrix non_hamiltonian_submatrix =
                        matrixBuilder.apply_to_subentity(space_.blocks[block],
                                                         pair.second);
                spectra_.at(pair.first).blocks[block] =
                        spectrumBuilder.apply_to_subentity_non_energy(non_hamiltonian_submatrix,
                                                                      unitary_transformation_matrix);
            }
        }
    }
}

const Matrix &runner::Runner::getMatrix(QuantityEnum quantity_enum) const {
    return matrices_.at(quantity_enum);
}

const Spectrum &runner::Runner::getSpectrum(QuantityEnum quantity_enum) const {
    return spectra_.at(quantity_enum);
}
