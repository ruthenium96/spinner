
#include <entities/Space.h>

#include "Space.h"

Space::Space(uint32_t total_space_size) {
    blocks.emplace_back();
    Subspace& lex_block = blocks[0];
    lex_block.basis.resize(total_space_size);
    for (uint32_t lex = 0; lex < total_space_size; ++lex) {
        lex_block.basis[lex][lex] = 1.0;
    }
}

Space::Space(std::vector<Subspace>&& v, Space::History h) {

    for (auto& subspace : v) {
        if (!subspace.basis.empty()) {
            blocks.emplace_back(std::move(subspace));
        }
    }

    history = h;
}

std::ostream &operator<<(std::ostream &os, const Space &space) {
    for (const Subspace& subspace : space.blocks) {
        os << subspace;
    }
    os << "------" << std::endl;
    return os;
}
