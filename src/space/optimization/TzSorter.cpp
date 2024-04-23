#include "TzSorter.h"

#include <cassert>
#include <utility>

namespace space::optimization {

Space TzSorter::apply(Space&& space) const {
    assert(!space.getBlocks().empty());
    auto totalSpaceSize = space.getBlocks()[0].decomposition->size_rows();
    std::vector<Subspace> vector_result;
    vector_result.reserve(space.getBlocks().size() * (max_ntz_proj + 1));
    for (size_t i = 0; i < space.getBlocks().size() * (max_ntz_proj + 1); ++i) {
        vector_result.emplace_back(factories_.createSparseSemiunitaryMatrix(0, totalSpaceSize));
    }

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace_parent = space.getBlocks()[i];

        for (size_t ntz_proj = 0; ntz_proj < max_ntz_proj + 1; ++ntz_proj) {
            BlockProperties block_properties = subspace_parent.properties;
            block_properties.n_proj = ntz_proj;
            vector_result[(max_ntz_proj + 1) * i + ntz_proj].properties = block_properties;
        }

        for (uint32_t l = 0; l < subspace_parent.decomposition->size_cols(); ++l) {
            // Value of total projection is calculated from the first index of map.
            // NB: there is no validation of the fact, that all indexes of decomposition
            // correspond to the same projection value, user should check it yourself.
            uint32_t index = subspace_parent.decomposition->GetNewIterator(l)->getNext().index;
            uint8_t ntz_proj = converter_->convert_lex_index_to_tz_projection(index);
            size_t j = (max_ntz_proj + 1) * i + ntz_proj;
            vector_result[j].decomposition->move_vector_from(l, subspace_parent.decomposition);
        }

        subspace_parent.decomposition->clear();
    }

    return Space(std::move(vector_result));
}

TzSorter::TzSorter(
    std::shared_ptr<const lexicographic::IndexConverter> indexes,
    quantum::linear_algebra::FactoriesList factories) :
    converter_(std::move(indexes)),
    factories_(std::move(factories)) {
    max_ntz_proj = converter_->get_max_ntz_proj();
}
}  // namespace space::optimization