#include "Subspectrum.h"

#include <utility>

std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum) {
    os << subspectrum.properties << std::endl;
    subspectrum.raw_data->print(os);
    os << std::endl;
    return os;
}

std::pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractMatrix>>
Subspectrum::energy(const Submatrix& hamiltonian_submatrix) {
    auto eigencouple = hamiltonian_submatrix.raw_data->diagonalizeValuesVectors();

    auto energy_subspectrum =
        Subspectrum(std::move(eigencouple.eigenvalues), hamiltonian_submatrix.properties);

    std::pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractMatrix>> answer = {
        std::move(energy_subspectrum),
        std::move(eigencouple.eigenvectors)};
    return answer;
}

Subspectrum Subspectrum::non_energy(
    const Submatrix& non_hamiltonian_submatrix,
    const std::unique_ptr<quantum::linear_algebra::AbstractMatrix>& unitary_transformation_matrix) {
    // TODO: I guess, we can do it faster. If B = U * A * U^T,
    //  we need only B_{ii} = \sum_{k} \sum_{l} U_{ik}*U_{il}*A_{kl}
    //  We can multiply A * U^T directly for O(N^3), then calculate all B_{ii} for O(N^2).
    auto raw_data =
        unitary_transformation_matrix->unitary_transform(non_hamiltonian_submatrix.raw_data)
            ->return_main_diagonal();

    auto non_energy_subspectrum =
        Subspectrum(std::move(raw_data), non_hamiltonian_submatrix.properties);

    return non_energy_subspectrum;
}
Subspectrum::Subspectrum(
    std::unique_ptr<quantum::linear_algebra::AbstractVector> raw_data_,
    BlockProperties properties_) {
    raw_data = std::move(raw_data_);
    properties = std::move(properties_);
}
