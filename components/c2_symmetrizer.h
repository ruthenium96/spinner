/* TODO: 
  * use #pragma once
  * use using for long type names
  * use int aliases like uint32 etc.
  * use proper naming 

*/
#ifndef JULY_C2_SYMMETRIZER_H
#define JULY_C2_SYMMETRIZER_H

#include "cmath"
#include "common/Task.h"
#include "unordered_map"
#include "unordered_set"
#include <boost/functional/hash.hpp>
#include <vector>

class Symmetrizer {
  public:
    Symmetrizer(std::vector<int> mults_, int pairs_);
    Task& operator()(Task& T);

    unsigned long symmetrized_lex(const unsigned long lex) const;

    std::vector<std::map<Index, Coefficient>>
    projector(std::map<Index, Coefficient>& m,
              std::unordered_map<size_t, size_t>& hs);

    static void add_to_hash_table(std::map<Index, Coefficient>& m,
                                  std::unordered_map<size_t, size_t>& hs);

    static void
    erase_if_zero(std::vector<std::map<Index, Coefficient>>& projections);

    static bool is_in_hash_table(const std::map<Index, Coefficient>& m,
                                 std::unordered_map<size_t, size_t>& hs);

  private:
    const std::vector<int> mults;
    std::vector<int> cumulative_product;
    int pairs;
    int tensor_size;
    int max_repr;
};

#endif // JULY_C2_SYMMETRIZER_H
