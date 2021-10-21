#ifndef JULY_MATRIXBUILDER_H
#define JULY_MATRIXBUILDER_H

#include <armadillo>
#include <common/LexicographicIndexConverter.h>
#include "entities/matrix/Matrix.h"
#include <entities/space/Space.h>
#include <entities/operator/Operator.h>

class MatrixBuilder {
public:
    explicit MatrixBuilder(spaces::LexicographicIndexConverter converter);

    Matrix apply(const Space& space, const Operator& new_operator);
private:
    const spaces::LexicographicIndexConverter converter_;
    Submatrix apply_to_subentity(const Subspace& subspace, const Operator& new_operator);
};
#endif //JULY_MATRIXBUILDER_H
