#ifndef JULY_INDEXES_H
#define JULY_INDEXES_H

#include <functional>
#include <numeric>
#include <vector>
#include <cstdint>

namespace spaces {
    class LexicographicIndexWorker {
    public:
        explicit LexicographicIndexWorker(std::vector<int> mults);
        uint8_t lex_to_ntzproj (uint32_t lex) const;
        std::vector<uint8_t> lex_to_nzs(uint32_t lex) const;
        uint32_t nzs_to_lex(const std::vector<uint8_t>& nzs) const;

        uint32_t total_space_size;
        const std::vector<int> mults_;
    private:
        std::vector<uint32_t> cumulative_product;
    };
}


#endif //JULY_INDEXES_H
