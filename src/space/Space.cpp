#include "Space.h"
namespace space {
Space::Space(uint32_t total_space_size) {
    UnitarySparseMatrix identity_decomposition;
    identity_decomposition.tensor_size = total_space_size;
    identity_decomposition.resize(total_space_size);
    for (uint32_t lex = 0; lex < total_space_size; ++lex) {
        identity_decomposition.add_to_position(1.0, lex, lex);
    }
    blocks_.emplace_back(std::move(identity_decomposition));
}

Space::Space(std::vector<Subspace>&& v) {
    for (auto& subspace : v) {
        if (!subspace.decomposition.empty()) {
            blocks_.emplace_back(std::move(subspace));
        }
    }
}

std::vector<Subspace>& Space::getBlocks() {
    return blocks_;
}

const std::vector<Subspace>& Space::getBlocks() const {
    return blocks_;
}

}  // namespace space

std::ostream& operator<<(std::ostream& os, const space::Space& space) {
    for (const space::Subspace& subspace : space.getBlocks()) {
        os << subspace;
    }
    os << "------" << std::endl;
    return os;
}
