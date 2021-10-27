#include "SpectrumBuilder.h"

Spectrum SpectrumBuilder::apply(const Matrix &hamiltonian_matrix, const std::vector<Matrix> &non_hamiltonian_matrices) {

    std::vector<Subspectrum> vector_result;
    vector_result.resize(hamiltonian_matrix.blocks.size());

    for (size_t i = 0; i < hamiltonian_matrix.blocks.size(); ++i) {
        DenseVector energy_values;
        DenseMatrix unitary_transformation_matrix;
        hamiltonian_matrix.blocks[i].raw_data.diagonalize(energy_values, unitary_transformation_matrix);
        vector_result[i].energy_raw_data = std::move(energy_values);
        for (const Matrix& non_hamiltonian_matrix : non_hamiltonian_matrices) {
            // TODO: I guess, we can do it faster. If B = U * A * U^T,
            //  we need only B_{ii} = \sum_{k} \sum_{l} U_{ik}*U_{il}*A_{kl}
            //  We can multiply A * U^T directly for O(N^3), then calculate all B_{ii} for O(N^2).
            DenseVector non_energy_values = unitary_transformation_matrix
                    .unitary_transform(non_hamiltonian_matrix.blocks[i].raw_data)
                    .return_main_diagonal();
            vector_result[i].non_energy_raw_data.emplace_back(
                    std::move(non_energy_values)
                    );
        }
    }
    return Spectrum(std::move(vector_result));
}
