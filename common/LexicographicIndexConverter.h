#ifndef JULY_INDEXES_H
#define JULY_INDEXES_H

#include <functional>
#include <numeric>
#include <vector>
#include <cstdint>

namespace spaces {
    class LexicographicIndexConverter {
    public:
        explicit LexicographicIndexConverter(std::vector<int> mults);
        uint8_t convert_lex_index_to_tz_projection (uint32_t lex) const;
        std::vector<uint8_t> convert_lex_index_to_sz_projections(uint32_t lex) const;
        uint32_t convert_sz_projections_to_lex_index(const std::vector<uint8_t>& nzs) const;

        uint32_t total_space_size;
        const std::vector<int> mults_;
    private:
        std::vector<uint32_t> cumulative_product;
    };
}


#endif //JULY_INDEXES_H
