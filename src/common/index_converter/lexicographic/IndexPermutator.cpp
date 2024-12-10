#include "IndexPermutator.h"
#include "src/group/Group.h"

namespace index_converter::lexicographic {

IndexPermutator::IndexPermutator(std::shared_ptr<const IndexConverter> converter, const group::Group& group) :
	converter_(converter), group_(group) {

}

uint32_t IndexPermutator::get_total_space_size() const {
	return converter_->get_total_space_size();
}

std::vector<IndexWithSign> IndexPermutator::convert_index_to_permutated_indexes(uint32_t index) const {
    std::vector<uint8_t> nzs = converter_->convert_lex_index_to_all_sz_projections(index);

    std::vector<std::vector<uint8_t>> permutated_vectors(group_.getElements().size());
    for (size_t i = 0; i < group_.getElements().size(); ++i) {
        permutated_vectors[i].resize(nzs.size());
        for (size_t j = 0; j < nzs.size(); ++j) {
            size_t position = group_.getElements()[i][j];
            permutated_vectors[i][position] = nzs[j];
        }
    }

    std::vector<IndexWithSign> answer;
    answer.resize(permutated_vectors.size());

    for (uint8_t g = 0; g < permutated_vectors.size(); ++g) {
        answer[g].index = converter_->convert_sz_projections_to_lex_index(permutated_vectors[g]);
        answer[g].sign = +1;
    }

    return answer;
}

} // namespace index_converter::lexicographic