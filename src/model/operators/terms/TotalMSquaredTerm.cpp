#include "TotalMSquaredTerm.h"

namespace model::operators {

TotalMSquaredTerm::TotalMSquaredTerm(
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter) :
    converter_(converter) {}

std::unique_ptr<Term> TotalMSquaredTerm::clone() const {
    return std::make_unique<TotalMSquaredTerm>(converter_);
};

void TotalMSquaredTerm::construct(
    quantum::linear_algebra::AbstractSymmetricMatrix&
        matrix_in_lexicografical_basis,
    const std::set<unsigned int>& indexes_of_vectors) const {

    for (const auto index_of_vector : indexes_of_vectors) {
        double max_spin = ((double)converter_->get_max_ntz_proj() - 1.0) / 2.0;
        double total_m = (double)converter_->convert_index_to_tz_projection(index_of_vector) 
            - max_spin;
        double total_m_squared = total_m * total_m;
        matrix_in_lexicografical_basis.add_to_position(
            total_m_squared,
            index_of_vector,
            index_of_vector);
    }
}

} // namespace model::operators