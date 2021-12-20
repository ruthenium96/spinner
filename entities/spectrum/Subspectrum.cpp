#include "Subspectrum.h"

std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum) {
    os << subspectrum.properties << std::endl;
    os << subspectrum.raw_data << std::endl;
    return os;
}

Subspectrum Subspectrum::energy(
    const Submatrix& hamiltonian_submatrix,
    DenseMatrix& unitary_transformation_matrix) {
    Subspectrum energy_subspectrum;
    hamiltonian_submatrix.raw_data.diagonalize(
        energy_subspectrum.raw_data,
        unitary_transformation_matrix);
    energy_subspectrum.properties = hamiltonian_submatrix.properties;
    return energy_subspectrum;
}

Subspectrum Subspectrum::non_energy(
    const Submatrix& non_hamiltonian_submatrix,
    const DenseMatrix& unitary_transformation_matrix) {
    Subspectrum non_energy_subspectrum;
    // TODO: I guess, we can do it faster. If B = U * A * U^T,
    //  we need only B_{ii} = \sum_{k} \sum_{l} U_{ik}*U_{il}*A_{kl}
    //  We can multiply A * U^T directly for O(N^3), then calculate all B_{ii} for O(N^2).
    non_energy_subspectrum.raw_data =
        unitary_transformation_matrix.unitary_transform(non_hamiltonian_submatrix.raw_data)
            .return_main_diagonal();
    non_energy_subspectrum.properties = non_hamiltonian_submatrix.properties;
    return non_energy_subspectrum;
}
