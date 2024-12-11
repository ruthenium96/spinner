#ifndef SPINNER_ABSTRACTINDEXPERMUTATOR_H
#define SPINNER_ABSTRACTINDEXPERMUTATOR_H

#include <cstdint>
#include <vector>


namespace index_converter {

struct IndexWithSign {
    uint32_t index;
    int8_t sign;
};

class AbstractIndexPermutator {
public:
    virtual std::vector<IndexWithSign> convert_index_to_permutated_indexes(uint32_t index) const = 0;
    virtual uint32_t get_total_space_size() const = 0;
    virtual ~AbstractIndexPermutator() = default;

};
} // namespace index_converter

#endif // SPINNER_ABSTRACTINDEXPERMUTATOR_H