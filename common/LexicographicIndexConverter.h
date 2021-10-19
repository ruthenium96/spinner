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
        [[nodiscard]] uint8_t convert_lex_index_to_tz_projection (uint32_t lex) const;
        [[nodiscard]] std::vector<uint8_t> convert_lex_index_to_all_sz_projections(uint32_t lex) const;
        [[nodiscard]] uint8_t convert_lex_index_to_one_sz_projection(uint32_t lex, uint32_t center) const;
        [[nodiscard]] uint32_t convert_sz_projections_to_lex_index(const std::vector<uint8_t>& nzs) const;

        // NB: one should check projection value before ladding
        [[nodiscard]] uint32_t ladder_projection(uint32_t lex, uint32_t center, int ladder) const;

        uint32_t total_space_size;
        [[nodiscard]] const std::vector<int>& get_mults() const;
        [[nodiscard]] const std::vector<double>& get_spins() const;

    private:
        std::vector<uint32_t> cumulative_product;
        const std::vector<int> mults_;
        std::vector<double> spins_;
    };
}


#endif //JULY_INDEXES_H
