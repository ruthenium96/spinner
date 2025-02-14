#include "TSquaredSorter.h"

#include <cassert>

namespace space::optimization {
TSquaredSorter::TSquaredSorter(std::shared_ptr<const index_converter::s_squared::IndexConverter> converter,
    quantum::linear_algebra::FactoriesList factories) :
	converter_(std::move(converter)),
    factories_(std::move(factories)) {
	max_total_mult_ = converter_->get_max_ntz_proj();
}

Space TSquaredSorter::apply(Space&& space) const {
    assert(!space.getBlocks().empty());
    auto totalSpaceSize = space.getBlocks()[0].decomposition->size_rows();

	size_t max_number_of_block = multiplicity_to_number_of_block(max_total_mult_);

    std::vector<Subspace> vector_result;
    vector_result.reserve(space.getBlocks().size() * (max_number_of_block + 1));
    for (size_t i = 0; i < space.getBlocks().size() * (max_number_of_block + 1); ++i) {
        vector_result.emplace_back(factories_.createSparseSemiunitaryMatrix(0, totalSpaceSize));
    }

#pragma omp parallel for shared(space, vector_result, max_number_of_block) default(none)
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace_parent = space.getBlocks()[i];

        for (size_t number_of_block = 0; number_of_block <= max_number_of_block; ++number_of_block) {
			spin_algebra::Multiplicity mult_of_block = number_of_block_to_multiplicity(number_of_block);
            BlockProperties block_properties = subspace_parent.properties;
            block_properties.total_mult = mult_of_block;
            vector_result[(max_number_of_block + 1) * i + number_of_block].properties = block_properties;
        }

        for (uint32_t l = 0; l < subspace_parent.decomposition->size_cols(); ++l) {
            // Value of total projection is calculated from the first index of map.
            // NB: there is no validation of the fact, that all indexes of decomposition
            // correspond to the same projection value, user should check it yourself.
            uint32_t index = subspace_parent.decomposition->GetNewIterator(l)->getNext().index;
            auto multiplicity = converter_->convert_index_to_total_multiplicity(index);

            size_t number_of_block = multiplicity_to_number_of_block(multiplicity);

            size_t j = (max_number_of_block + 1) * i + number_of_block;
            vector_result[j].decomposition->move_vector_from(l, subspace_parent.decomposition);
        }

        subspace_parent.decomposition->clear();
    }

    return Space(std::move(vector_result));
}


spin_algebra::Multiplicity TSquaredSorter::number_of_block_to_multiplicity(size_t number_of_block) const {
    if (max_total_mult_ % 2 == 0) {
        return 2 * number_of_block + 2;
    } else {
        return 2 * number_of_block + 1;
    }
}

size_t TSquaredSorter::multiplicity_to_number_of_block(spin_algebra::Multiplicity multiplicity) const {
    if (max_total_mult_ % 2 == 0) {
        return multiplicity / 2 - 1;
    } else {
        return (multiplicity - 1) / 2;
    }
}

}  // namespace space::optimization
