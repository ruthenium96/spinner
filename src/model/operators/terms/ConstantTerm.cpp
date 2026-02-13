#include "ConstantTerm.h"

#include <utility>

namespace spinner::model::operators {
void ConstantTerm::construct(
    linalg_structures::AbstractSymmetricMatrix&
        matrix_in_lexicografical_basis,
    const std::set<unsigned int>& indexes_of_vectors) const {
    for (const auto& index_of_vector : indexes_of_vectors) {
        matrix_in_lexicografical_basis.add_to_position(*constant_, index_of_vector, index_of_vector);
    }
}

ConstantTerm::ConstantTerm(std::shared_ptr<const double> constant) :
    constant_(std::move(constant)) {}

std::unique_ptr<Term> ConstantTerm::clone() const {
    return std::make_unique<ConstantTerm>(constant_);
}
} // namespace spinner::model::operators