#include "Spectrum.h"

std::ostream& operator<<(std::ostream& os, const Spectrum& spectrum) {
    for (const Subspectrum& subspectrum : spectrum.blocks) {
        os << subspectrum;
    }
    os << "------" << std::endl;
    return os;
}

Spectrum::Spectrum(std::vector<Subspectrum>&& m) : blocks(std::move(m)) {}

Spectrum Spectrum::non_energy(
    const Matrix& non_hamiltonian_matrix,
    const std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>&
        unitary_transformation_matrices) {
    std::vector<Subspectrum> vector_result;
    vector_result.reserve(non_hamiltonian_matrix.blocks.size());

    if (non_hamiltonian_matrix.blocks.size() != unitary_transformation_matrices.size()) {
        throw std::length_error(
            "non_hamiltonian_matrix.blocks.size() != unitary_transformation_matrices.size()");
    }

    for (size_t i = 0; i < non_hamiltonian_matrix.blocks.size(); ++i) {
        vector_result.emplace_back(Subspectrum::non_energy(
            non_hamiltonian_matrix.blocks[i],
            unitary_transformation_matrices[i]));
    }
    return Spectrum(std::move(vector_result));
}

std::pair<
    Spectrum,
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
Spectrum::energy(const Matrix& hamiltonian_matrix) {
    std::vector<Subspectrum> vector_subspectrum_result;
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        vector_unitary_transformation_matrices_result;
    vector_subspectrum_result.reserve(hamiltonian_matrix.blocks.size());
    vector_unitary_transformation_matrices_result.reserve(hamiltonian_matrix.blocks.size());

    for (const auto& hamiltonian_submatrix : hamiltonian_matrix.blocks) {
        auto [subspectrum_energy, unitary_transformation_matrix] =
            Subspectrum::energy(hamiltonian_submatrix);
        vector_subspectrum_result.emplace_back(std::move(subspectrum_energy));
        vector_unitary_transformation_matrices_result.emplace_back(
            std::move(unitary_transformation_matrix));
    }

    std::pair<
        Spectrum,
        std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
        answer = {
            Spectrum(std::move(vector_subspectrum_result)),
            std::move(vector_unitary_transformation_matrices_result)};

    return answer;
}
