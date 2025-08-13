#ifndef SPINNER_S2INDEXPERMUTATOR_H
#define SPINNER_S2INDEXPERMUTATOR_H

#include "src/common/index_converter/AbstractIndexPermutator.h"
#include "IndexConverter.h"
#include "src/group/Group.h"

namespace index_converter::s_squared {

class IndexPermutator : public AbstractIndexPermutator {
  public:
    IndexPermutator(std::shared_ptr<const IndexConverter> converter, const group::Group& group);
    std::vector<IndexWithSign> convert_index_to_permutated_indexes(uint32_t index) const override;
    uint32_t get_total_space_size() const override;

  private:
    std::shared_ptr<const IndexConverter> converter_;
    group::Permutation extendPermutation(const group::Permutation& permutation) const;
    std::vector<group::Permutation> extended_permutations_;
    size_t initial_size_of_permutation_;
    void construct_level_and_sign(const group::Permutation& extended_g, const Level& level, Level& permutated_level, int8_t& sign) const;
};
}

#endif // SPINNER_S2INDEXPERMUTATOR_H