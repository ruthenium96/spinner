#ifndef JULY_INDEXES_H
#define JULY_INDEXES_H

#include <functional>
#include <numeric>
#include <vector>

namespace Spaces {
    class Indexes {
    public:
        explicit Indexes(std::vector<int> mults);
        int lex_to_ntzproj (unsigned long lex) const;
        std::vector<int> lex_to_nzs(unsigned long lex) const;
        unsigned long nzs_to_lex(const std::vector<int>& nzs) const;

        unsigned int total_space_size;
        const std::vector<int> mults_;
    private:
        std::vector<int> cumulative_product;
    };
}


#endif //JULY_INDEXES_H
