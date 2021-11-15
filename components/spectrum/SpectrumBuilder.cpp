#include "SpectrumBuilder.h"

Spectrum SpectrumBuilder::apply_to_energy(const Matrix &hamiltonian_matrix) {
    std::vector<Subspectrum> vector_result;
    vector_result.resize(hamiltonian_matrix.blocks.size());
    unitary_transformation_matrices.resize(hamiltonian_matrix.blocks.size());

    for (size_t i = 0; i < hamiltonian_matrix.blocks.size(); ++i) {
        vector_result[i] = apply_to_subentity_energy(hamiltonian_matrix.blocks[i],
                                                     unitary_transformation_matrices[i]);
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
        vector_result[i] = apply_to_subentity_non_energy(non_hamiltonian_matrix.blocks[i],
                                                         unitary_transformation_matrices[i]);
    }
    return Spectrum(std::move(vector_result));
}

Subspectrum SpectrumBuilder::apply_to_subentity_energy(const Submatrix &hamiltonian_submatrix, DenseMatrix& unitary_transformation_matrix) {
    Subspectrum energy_subspectrum;
    hamiltonian_submatrix.raw_data.diagonalize(energy_subspectrum.raw_data, unitary_transformation_matrix);
    energy_subspectrum.properties = hamiltonian_submatrix.properties;
    return energy_subspectrum;
}

Subspectrum SpectrumBuilder::apply_to_subentity_non_energy(const Submatrix &non_hamiltonian_submatrix,
                                                           const DenseMatrix &unitary_transformation_matrix) {
    Subspectrum non_energy_subspectrum;
    // TODO: I guess, we can do it faster. If B = U * A * U^T,
    //  we need only B_{ii} = \sum_{k} \sum_{l} U_{ik}*U_{il}*A_{kl}
    //  We can multiply A * U^T directly for O(N^3), then calculate all B_{ii} for O(N^2).
    non_energy_subspectrum.raw_data = unitary_transformation_matrix
            .unitary_transform(non_hamiltonian_submatrix.raw_data).return_main_diagonal();
    non_energy_subspectrum.properties = non_hamiltonian_submatrix.properties;
    return non_energy_subspectrum;
}
