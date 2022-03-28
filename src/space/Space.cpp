#include "Space.h"
namespace space {
Space::Space(uint32_t total_space_size) {
    UnitarySparseMatrix identity_decomposition;
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
}  // namespace space

std::ostream& operator<<(std::ostream& os, const space::Space& space) {
    for (const space::Subspace& subspace : space.blocks) {
        os << subspace;
    }
    os << "------" << std::endl;
    return os;
}
