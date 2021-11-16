#ifndef JULY_MATRIXBUILDER_H
#define JULY_MATRIXBUILDER_H

#include <armadillo>
#include <common/lexicographic/IndexConverter.h>
#include "entities/matrix/Matrix.h"
#include <entities/space/Space.h>
#include <entities/operator/Operator.h>

class MatrixBuilder {
public:
    explicit MatrixBuilder(lexicographic::IndexConverter converter);

    Matrix apply(const Space& space, const Operator& new_operator);
    Submatrix apply_to_subentity(const Subspace& subspace, const Operator& new_operator);
private:
    const lexicographic::IndexConverter converter_;
};
#endif //JULY_MATRIXBUILDER_H
