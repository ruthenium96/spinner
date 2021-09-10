#include "LexicographicIndexConverter.h"

spaces::LexicographicIndexConverter::LexicographicIndexConverter(std::vector<int> mults) : mults_(mults) {
    // auxiliary vector for uint32_t <-> std::vector<uint8_t> transformation
    cumulative_product.resize(mults.size() + 1);
    cumulative_product.back() = 1;
    std::partial_sum(mults.rbegin(), mults.rend(),
                     cumulative_product.rbegin() + 1, std::multiplies<>());

    // It is the size of our projection-based space.
    total_space_size = cumulative_product[0];
}

uint8_t spaces::LexicographicIndexConverter::convert_lex_index_to_tz_projection(uint32_t lex) const {
    uint8_t ntz_proj = 0;
    for (int i = 0; i < mults_.size(); ++i) {
        ntz_proj += (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return ntz_proj;
}

std::vector<uint8_t> spaces::LexicographicIndexConverter::convert_lex_index_to_sz_projections(uint32_t lex) const {
    std::vector<uint8_t> nzs(mults_.size());
    for (int i = 0; i < mults_.size(); ++i) {
        nzs[i] = (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return(std::move(nzs));
}

uint32_t spaces::LexicographicIndexConverter::convert_sz_projections_to_lex_index(const std::vector<uint8_t> &nzs) const {
    uint32_t lex = 0;
    for (int i = 0; i < mults_.size(); ++i) {
        lex += nzs[i] * cumulative_product[i + 1];
    }
    return lex;
}
