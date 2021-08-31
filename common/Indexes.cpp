#include "Indexes.h"

Spaces::Indexes::Indexes(std::vector<int> mults) : mults_(mults) {
    // auxiliary vector for Lex_Index <-> std::vector<Projection> transformation
    cumulative_product.resize(mults.size() + 1);
    cumulative_product.back() = 1;
    std::partial_sum(mults.rbegin(), mults.rend(),
                     cumulative_product.rbegin() + 1, std::multiplies<>());

    // It is the size of our projection-based space.
    total_space_size = cumulative_product[0];
}

Projection Spaces::Indexes::lex_to_ntzproj(Lex_Index lex) const {
    Projection ntz_proj = 0;
    for (int i = 0; i < mults_.size(); ++i) {
        ntz_proj += (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return ntz_proj;
}

std::vector<Projection> Spaces::Indexes::lex_to_nzs(Lex_Index lex) const {
    std::vector<Projection> nzs(mults_.size());
    for (int i = 0; i < mults_.size(); ++i) {
        nzs[i] = (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return(std::move(nzs));
}

Lex_Index Spaces::Indexes::nzs_to_lex(const std::vector<Projection> &nzs) const {
    Lex_Index lex = 0;
    for (int i = 0; i < mults_.size(); ++i) {
        lex += nzs[i] * cumulative_product[i + 1];
    }
    return lex;
}
