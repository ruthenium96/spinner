#include "IndexConverter.h"

lexicographic::IndexConverter::IndexConverter(std::vector<int> mults) : mults_(mults) {
    // auxiliary vector for uint32_t <-> std::vector<uint8_t> transformation
    cumulative_product.resize(mults.size() + 1);
    cumulative_product.back() = 1;
    std::partial_sum(
        mults.rbegin(),
        mults.rend(),
        cumulative_product.rbegin() + 1,
        std::multiplies<>());

    // It is the size of our projection-based space.
    total_space_size = cumulative_product[0];

    // mult = 2 * spin + 1
    spins_.resize(mults_.size());
    for (size_t i = 0; i < mults.size(); ++i) {
        spins_[i] = (mults_[i] - 1) / 2.0;
    }
}

uint8_t lexicographic::IndexConverter::convert_lex_index_to_tz_projection(uint32_t lex) const {
    uint8_t ntz_proj = 0;
    for (size_t i = 0; i < mults_.size(); ++i) {
        ntz_proj += (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return ntz_proj;
}

std::vector<uint8_t>
lexicographic::IndexConverter::convert_lex_index_to_all_sz_projections(uint32_t lex) const {
    std::vector<uint8_t> nzs(mults_.size());
    for (size_t i = 0; i < mults_.size(); ++i) {
        nzs[i] = (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return nzs;
}

uint8_t lexicographic::IndexConverter::convert_lex_index_to_one_sz_projection(
    uint32_t lex,
    uint32_t center) const {
    return (lex % cumulative_product[center]) / cumulative_product[center + 1];
}

uint32_t lexicographic::IndexConverter::convert_sz_projections_to_lex_index(
    const std::vector<uint8_t>& nzs) const {
    uint32_t lex = 0;
    for (size_t i = 0; i < mults_.size(); ++i) {
        lex += nzs[i] * cumulative_product[i + 1];
    }
    return lex;
}

uint32_t
lexicographic::IndexConverter::ladder_projection(uint32_t lex, uint32_t center, int ladder) const {
    return lex + ladder * cumulative_product[center + 1];
}

const std::vector<int>& lexicographic::IndexConverter::get_mults() const {
    return mults_;
}

const std::vector<double>& lexicographic::IndexConverter::get_spins() const {
    return spins_;
}
uint32_t lexicographic::IndexConverter::get_total_space_size() const {
    return total_space_size;
}
