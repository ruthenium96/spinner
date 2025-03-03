#include "TotalSSquaredTerm.h"

namespace model::operators::ito {

TotalSSquaredTerm::TotalSSquaredTerm(
	std::shared_ptr<const index_converter::s_squared::IndexConverter> converter) :
	converter_(converter) {}

std::unique_ptr<Term> TotalSSquaredTerm::clone() const {
	return std::make_unique<TotalSSquaredTerm>(converter_);
}

void TotalSSquaredTerm::construct(
	quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
	const std::set<unsigned int>& indexes_of_vectors) const {
    for (const auto& index : indexes_of_vectors) {
		auto multiplicity = converter_->convert_index_to_total_multiplicity(index);
		auto spin = ((double)multiplicity - 1.0) / 2.0;
		matrix.add_to_position(spin * (spin + 1), index, index);
	}
}


} // namespace model::operators::ito