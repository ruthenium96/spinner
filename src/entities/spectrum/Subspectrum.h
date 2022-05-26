#ifndef SPINNER_SUBSPECTRUM_H
#define SPINNER_SUBSPECTRUM_H

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/DenseMatrix.h"
#include "src/entities/matrix/Submatrix.h"

struct Subspectrum {
    BlockProperties properties;
    DenseVector raw_data;

    static Subspectrum
    energy(const Submatrix& hamiltonian_submatrix, DenseMatrix& unitary_transformation_matrix);
    static Subspectrum non_energy(
        const Submatrix& non_hamiltonian_submatrix,
        const DenseMatrix& unitary_transformation_matrix);

    friend std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum);
};

#endif  //SPINNER_SUBSPECTRUM_H
