#ifndef SPINNER_BLOCKPROPERTIES_H
#define SPINNER_BLOCKPROPERTIES_H

#include <optional>
#include <ostream>
#include <string>
#include <vector>

struct BlockProperties {
    std::string get_representation_name() const;

    std::optional<uint32_t> n_proj = std::nullopt;
    uint32_t dimensionality = 1;
    uint32_t degeneracy = 1;
    std::vector<uint32_t> representation;

    bool operator==(const BlockProperties&) const = default;
};

// TODO: consider moving all the printing code out of the domain classes
std::ostream& operator<<(std::ostream& os, const BlockProperties& properties);

#endif  //SPINNER_BLOCKPROPERTIES_H
