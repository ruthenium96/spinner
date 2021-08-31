#ifndef JULY_INDEXES_H
#define JULY_INDEXES_H

#include <functional>
#include <numeric>
#include <vector>
#include "Subspace.h"

using Projection = uint8_t;

namespace Spaces {
    class Indexes {
    public:
        explicit Indexes(std::vector<int> mults);
        Projection lex_to_ntzproj (Lex_Index lex) const;
        std::vector<Projection> lex_to_nzs(Lex_Index lex) const;
        Lex_Index nzs_to_lex(const std::vector<Projection>& nzs) const;

        uint32_t total_space_size;
        const std::vector<int> mults_;
    private:
        std::vector<uint32_t> cumulative_product;
    };
}


#endif //JULY_INDEXES_H
