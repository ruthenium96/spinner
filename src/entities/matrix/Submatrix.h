#ifndef SPINNER_SUBMATRIX_H
#define SPINNER_SUBMATRIX_H

#include <src/entities/data_structures/FactoriesList.h>

#include "src/entities/BlockProperties.h"
#include "src/model/operators/Operator.h"
#include "src/space/Subspace.h"

#include "src/common/index_converter/AbstractIndexConverter.h"

namespace spinner {

struct Submatrix {
    BlockProperties properties;
    std::unique_ptr<linear_algebra::AbstractDiagonalizableMatrix> raw_data;

    Submatrix() = default;

    Submatrix(
        std::unique_ptr<linear_algebra::AbstractDiagonalizableMatrix> raw_data_,
        BlockProperties properties_);

    Submatrix(
        const space::Subspace& subspace,
        const model::operators::Operator& new_operator,
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        const linear_algebra::FactoriesList& factories,
        bool return_sparse_if_possible);
};

} // namespace spinner

#endif  //SPINNER_SUBMATRIX_H
