#include "SpectrumBuilder.h"

Spectrum SpectrumBuilder::apply_to_energy(const Matrix &hamiltonian_matrix) {
    std::vector<Subspectrum> vector_result;
    vector_result.resize(hamiltonian_matrix.blocks.size());
    unitary_transformation_matrices.resize(hamiltonian_matrix.blocks.size());

    for (size_t i = 0; i < hamiltonian_matrix.blocks.size(); ++i) {
        DenseVector energy_values;
        DenseMatrix& unitary_transformation_matrix = unitary_transformation_matrices[i];
        hamiltonian_matrix.blocks[i].raw_data.diagonalize(energy_values, unitary_transformation_matrix);
        vector_result[i].raw_data = std::move(energy_values);
        vector_result[i].properties = hamiltonian_matrix.blocks[i].properties;
    }
    return Spectrum(std::move(vector_result));
}

Spectrum SpectrumBuilder::apply_to_non_energy(const Matrix &non_hamiltonian_matrix) {
    std::vector<Subspectrum> vector_result;
    vector_result.resize(non_hamiltonian_matrix.blocks.size());

    if (vector_result.size() != unitary_transformation_matrices.size()) {
        throw std::length_error("vector_result.size() != unitary_transformation_matrices.size()");
    }

    for (size_t i = 0; i < non_hamiltonian_matrix.blocks.size(); ++i) {
        // TODO: I guess, we can do it faster. If B = U * A * U^T,
        //  we need only B_{ii} = \sum_{k} \sum_{l} U_{ik}*U_{il}*A_{kl}
        //  We can multiply A * U^T directly for O(N^3), then calculate all B_{ii} for O(N^2).
        DenseVector non_energy_values = unitary_transformation_matrices[i]
                .unitary_transform(non_hamiltonian_matrix.blocks[i].raw_data)
                .return_main_diagonal();
                vector_result[i].raw_data = std::move(non_energy_values);
                vector_result[i].properties = non_hamiltonian_matrix.blocks[i].properties;
    }
    return Spectrum(std::move(vector_result));
}
