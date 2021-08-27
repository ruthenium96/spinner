/* TODO: 
  * use using for long type names
  * use int aliases like uint32 etc.
*/
#ifndef JULY_C2_SYMMETRIZER_H
#define JULY_C2_SYMMETRIZER_H

#include "cmath"
#include "common/Space.h"
#include "unordered_map"
#include <boost/functional/hash.hpp>
#include <vector>

class Symmetrizer {
  public:
    Symmetrizer(const Spaces::Indexes& indexes, int pairs_);
    Space& operator()(Space& space);

    unsigned long symmetrized_lex(unsigned long lex) const;

    std::vector<Decomposition> projector(Decomposition& m, std::unordered_map<size_t, size_t>& hs);

    static void add_to_hash_table(Decomposition& m, std::unordered_map<size_t, size_t>& hs);

    static void erase_if_zero(std::vector<Decomposition>& projections);

    static bool is_in_hash_table(const Decomposition& m, std::unordered_map<size_t, size_t>& hs);

  private:
    const Spaces::Indexes& indexes_;
    int pairs;
    int max_repr;
};

#endif // JULY_C2_SYMMETRIZER_H
