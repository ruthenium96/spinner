#include "ScalarProduct.h"

#include <cmath>

ScalarProduct::ScalarProduct(size_t number_of_spins) {
    std::shared_ptr<DenseMatrix> mutable_coefficients = std::make_unique<DenseMatrix>();
    mutable_coefficients->resize_with_nans(number_of_spins, number_of_spins);
    for (size_t i = 0; i < number_of_spins; ++i) {
        for (size_t j = i + 1; j < number_of_spins; ++j) {
            // TODO: it does not looks good
            mutable_coefficients->assign_to_position(-1, i, j);
            mutable_coefficients->assign_to_position(-1, j, i);
        }
    }
    coefficients = mutable_coefficients;
}

ScalarProduct::ScalarProduct(std::shared_ptr<const DenseMatrix> parameters) :
    coefficients(std::move(parameters)) {}

void ScalarProduct::construct(
    lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b) const {
    if (!std::isnan(coefficients->operator()(center_a, center_b))) {
        double factor = -2 * coefficients->operator()(center_a, center_b);
        matrix_in_lexicografical_basis
            .add_scalar_product(index_of_vector, center_a, center_b, factor);
    }
}

std::shared_ptr<const DenseMatrix> ScalarProduct::get_parameters() const {
    return coefficients;
}
