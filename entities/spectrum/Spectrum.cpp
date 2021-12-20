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
    const std::vector<DenseMatrix>& unitary_transformation_matrices) {
    std::vector<Subspectrum> vector_result;
    vector_result.resize(non_hamiltonian_matrix.blocks.size());

    if (vector_result.size() != unitary_transformation_matrices.size()) {
        throw std::length_error("vector_result.size() != unitary_transformation_matrices.size()");
    }

    for (size_t i = 0; i < non_hamiltonian_matrix.blocks.size(); ++i) {
        vector_result[i] = Subspectrum::non_energy(
            non_hamiltonian_matrix.blocks[i],
            unitary_transformation_matrices[i]);
    }
    return Spectrum(std::move(vector_result));
}

Spectrum Spectrum::energy(
    const Matrix& hamiltonian_matrix,
    std::vector<DenseMatrix>& unitary_transformation_matrices) {
    std::vector<Subspectrum> vector_result;
    vector_result.resize(hamiltonian_matrix.blocks.size());
    unitary_transformation_matrices.resize(hamiltonian_matrix.blocks.size());

    for (size_t i = 0; i < hamiltonian_matrix.blocks.size(); ++i) {
        vector_result[i] =
            Subspectrum::energy(hamiltonian_matrix.blocks[i], unitary_transformation_matrices[i]);
    }
    return Spectrum(std::move(vector_result));
}
