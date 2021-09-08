
#include <entities/Space.h>

#include "Space.h"

using namespace entities;

Space::Space(uint32_t total_space_size): Entity(Entity::SPACE) {
    blocks.emplace_back();
    Subspace& lex_block = blocks[0];
    lex_block.basis.resize(total_space_size);
    for (uint32_t lex = 0; lex < total_space_size; ++lex) {
        lex_block.basis[lex][lex] = 1.0;
    }
}

std::ostream &operator<<(std::ostream &os, const Space &space) {
    for (const Subspace& Ss : space.blocks) {
        for (auto& m: Ss.basis) {
            for (auto& d: m) {
                os << d.second << "*[" << d.first << "] ";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << "------" << std::endl;
    return os;
}

Space::Space(std::vector<Subspace> v, entities::Entity::History h): Entity(Entity::SPACE) {

    blocks = v;
    history = h;
}
