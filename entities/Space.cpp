#include "Space.h"

Space::Space(uint32_t total_space_size) {
    Subspace lex_block;
    // TODO: do it better
    lex_block.tensor_size = total_space_size;
    lex_block.resize(total_space_size);
    for (uint32_t lex = 0; lex < total_space_size; ++lex) {
//        lex_block(lex, lex) = 1.0;
        lex_block.add_to_position(1.0, lex, lex);
    }
    blocks.emplace_back(std::move(lex_block));
}

Space::Space(std::vector<Subspace>&& v) {

    for (auto& subspace : v) {
        if (!subspace.empty()) {
            blocks.emplace_back(std::move(subspace));
        }
    }
}

std::ostream &operator<<(std::ostream &os, const Space &space) {
    for (const Subspace& subspace : space.blocks) {
        os << subspace;
    }
    os << "------" << std::endl;
    return os;
}
