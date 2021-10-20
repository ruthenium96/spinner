#include "Space.h"

Space::Space(uint32_t total_space_size) {
    NewBasisDecomposition identity_decomposition;
    identity_decomposition.tensor_size = total_space_size;
    identity_decomposition.resize(total_space_size);
    for (uint32_t lex = 0; lex < total_space_size; ++lex) {
        identity_decomposition.add_to_position(1.0, lex, lex);
    }
    blocks.emplace_back(std::move(identity_decomposition));
}

Space::Space(std::vector<Subspace>&& v) {

    for (auto& subspace : v) {
        if (!subspace.decomposition.empty()) {
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