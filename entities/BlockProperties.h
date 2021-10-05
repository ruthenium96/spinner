#ifndef JULY_BLOCKPROPERTIES_H
#define JULY_BLOCKPROPERTIES_H

#include <optional>
#include <ostream>
#include <string>
#include <vector>

struct BlockProperties {
    [[nodiscard]] std::string get_representation_name() const;

    std::optional<uint32_t> n_proj = std::nullopt;
    uint32_t dimensionality = 1;
    uint32_t degeneracy = 1;
    std::vector<uint32_t> representation;
};

std::ostream &operator<<(std::ostream &os, const BlockProperties &properties);

#endif //JULY_BLOCKPROPERTIES_H
