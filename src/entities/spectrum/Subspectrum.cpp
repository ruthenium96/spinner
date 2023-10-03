#include "Subspectrum.h"

#include <utility>

std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum) {
    os << subspectrum.properties << std::endl;
    subspectrum.raw_data->print(os);
    os << std::endl;
    return os;
}

std::pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
Subspectrum::energy(const Submatrix& hamiltonian_submatrix) {
    auto eigencouple = hamiltonian_submatrix.raw_data->diagonalizeValuesVectors();

    auto energy_subspectrum =
        Subspectrum(std::move(eigencouple.eigenvalues), hamiltonian_submatrix.properties);

    std::pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        answer = {std::move(energy_subspectrum), std::move(eigencouple.eigenvectors)};
    return answer;
}

Subspectrum Subspectrum::non_energy(
    const Submatrix& non_hamiltonian_submatrix,
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>&
        unitary_transformation_matrix) {
    auto raw_data = unitary_transformation_matrix->unitaryTransformAndReturnMainDiagonal(
        non_hamiltonian_submatrix.raw_data);

    auto non_energy_subspectrum =
        Subspectrum(std::move(raw_data), non_hamiltonian_submatrix.properties);

    return non_energy_subspectrum;
}
Subspectrum::Subspectrum(
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> raw_data_,
    BlockProperties properties_) {
    raw_data = std::move(raw_data_);
    properties = std::move(properties_);
}

Subspectrum Subspectrum::energy_without_eigenvectors(const Submatrix& hamiltonian_submatrix) {
    auto eigenvalues = hamiltonian_submatrix.raw_data->diagonalizeValues();

    auto energy_subspectrum = Subspectrum(std::move(eigenvalues), hamiltonian_submatrix.properties);

    return std::move(energy_subspectrum);
}
