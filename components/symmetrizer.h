#ifndef JULY_SYMMETRIZER_H
#define JULY_SYMMETRIZER_H

#include "cmath"
#include "common/Space.h"
#include "unordered_map"
#include <boost/functional/hash.hpp>
#include <vector>
#include <groups/Group.h>

class Symmetrizer {
  public:
    Symmetrizer(const Spaces::Indexes& indexes, const Group& group);
    Space& operator()(Space& space) const;

    std::vector<Decomposition> get_symmetrical_projected_decompositions(Decomposition& m, std::unordered_map<size_t, size_t>& hs) const;

    // TODO: these functions are not about symmetrization.
    //  Should we refactor them and create a new class?
    static void add_to_hash_table(Decomposition& m, std::unordered_map<size_t, size_t>& hs);

    static void erase_if_zero(std::vector<Decomposition>& projections);

    static bool is_in_hash_table(const Decomposition& m, std::unordered_map<size_t, size_t>& hs);

  private:
    const Spaces::Indexes& indexes_;
    const Group& group_;
};

#endif // JULY_SYMMETRIZER_H
