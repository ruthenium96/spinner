#ifndef JULY_SUBSPECTRUM_H
#define JULY_SUBSPECTRUM_H

#include "entities/BlockProperties.h"
#include "entities/data_structures/DenseMatrix.h"
#include "entities/matrix/Submatrix.h"

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

#endif  //JULY_SUBSPECTRUM_H
