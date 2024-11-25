#include "IndexConverter.h"

namespace index_converter::lexicographic {

IndexConverter::IndexConverter(std::vector<spin_algebra::Multiplicity> mults) :
    AbstractIndexConverter(mults) {
    // auxiliary vector for uint32_t <-> std::vector<uint8_t> transformation
    cumulative_product_.resize(get_mults().size() + 1);
    cumulative_product_.back() = 1;
    std::partial_sum(
        get_mults().rbegin(),
        get_mults().rend(),
        cumulative_product_.rbegin() + 1,
        std::multiplies<>());
}

std::vector<uint8_t> IndexConverter::convert_lex_index_to_all_sz_projections(uint32_t lex) const {
    std::vector<uint8_t> nzs(get_mults().size());
    for (size_t i = 0; i < get_mults().size(); ++i) {
        nzs[i] = (lex % cumulative_product_[i]) / cumulative_product_[i + 1];
    }
    return nzs;
}

uint8_t IndexConverter::convert_lex_index_to_one_sz_projection(
    uint32_t lex,
    uint32_t center) const {
    return (lex % cumulative_product_[center]) / cumulative_product_[center + 1];
}

uint32_t IndexConverter::convert_sz_projections_to_lex_index(
    const std::vector<uint8_t>& nzs) const {
    uint32_t lex = 0;
    for (size_t i = 0; i < get_mults().size(); ++i) {
        lex += nzs[i] * cumulative_product_[i + 1];
    }
    return lex;
}

std::vector<uint32_t>
IndexConverter::convert_index_to_permutated_indexes(
    uint32_t index, const group::Group& group) const {
    std::vector<uint8_t> nzs = convert_lex_index_to_all_sz_projections(index);
    std::vector<std::vector<uint8_t>> permutated_vectors = group.permutate(nzs);

    std::vector<uint32_t> answer;
    answer.resize(permutated_vectors.size());

    for (uint8_t g = 0; g < permutated_vectors.size(); ++g) {
        answer[g] = convert_sz_projections_to_lex_index(permutated_vectors[g]);
    }

    return std::move(answer);
}

uint8_t IndexConverter::convert_index_to_tz_projection(
    uint32_t index) const {
    uint8_t ntz_proj = 0;
    for (size_t i = 0; i < get_mults().size(); ++i) {
        ntz_proj += (index % cumulative_product_[i]) / cumulative_product_[i + 1];
    }
    return ntz_proj;
}

uint32_t
IndexConverter::ladder_projection(uint32_t lex, uint32_t center, int ladder) const {
    return lex + ladder * cumulative_product_[center + 1];
}

}