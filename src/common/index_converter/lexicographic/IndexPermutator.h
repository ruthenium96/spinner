#ifndef SPINNER_LEXINDEXPERMUTATOR_H
#define SPINNER_LEXINDEXPERMUTATOR_H

#include <memory>

#include "IndexConverter.h"
#include "src/common/index_converter/AbstractIndexPermutator.h"
#include "src/group/Group.h"

namespace index_converter::lexicographic {

class IndexPermutator : public AbstractIndexPermutator {
public:
    IndexPermutator(std::shared_ptr<const IndexConverter> converter, const group::Group& group);
    std::vector<IndexWithSign> convert_index_to_permutated_indexes(uint32_t index) const override;
    uint32_t get_total_space_size() const override;

private:
    std::shared_ptr<const IndexConverter> converter_;
    group::Group group_;
};
} // namespace index_converter::lexicographic

#endif // SPINNER_LEXINDEXPERMUTATOR_H