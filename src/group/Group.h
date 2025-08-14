#ifndef SPINNER_GROUP_H
#define SPINNER_GROUP_H

#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace group {

using PermutationArrayType = uint8_t;
using Permutation = std::vector<PermutationArrayType>;
using CayleyTable = std::map<std::pair<uint8_t, uint8_t>, std::set<uint8_t>>;

class Group {
  public:
    /*
     group_size : the number of permutations in the group

     number_of_representations : also the number of conjugacy classes

     is_abelian: we say that group is Abelian if group_size equals number_of_representations

     number_of_generators : the minimal number of elements required to produce all elements of the group

     group_in_form_of_generators : all elements of the group in the form a^vector[0] * b^vector[1] * ...,
                                   where {a, b, ...} -- generators in the _descending_ order of element order.
                                   element order is the minimum m, such that g^m = e
                                   It has the size (group_size) x (number of generators).

     coefficients_of_projectors : has the size (number_of_representations) x (dimension of representation) x (group_size).
                                  important: coefficients_of_projectors[0][0] should correspond to full-symmetric representation.

     cayley_table_table : table of direct products of representations
     */

    struct AlgebraicProperties {
        uint8_t group_size;
        uint8_t number_of_representations;
        bool is_abelian;
        std::vector<uint8_t> dimension_of_representation;
        std::vector<uint8_t> number_of_projectors_of_representation;
        uint8_t number_of_generators;
        std::vector<size_t> orders_of_generators;
        std::vector<std::vector<size_t>> group_in_form_of_generators;
        std::vector<size_t> orders_of_elements;
        std::vector<std::vector<std::vector<double>>> coefficients_of_projectors;
        CayleyTable cayley_table;
    };

    enum GroupTypeEnum {
        S2,
        Dihedral,
    };

    struct GroupType {
        GroupTypeEnum type_enum;
        std::optional<unsigned int> order;

        GroupType(GroupTypeEnum enum_) : type_enum(enum_), order(std::nullopt) {};
        GroupType(GroupTypeEnum enum_, std::optional<unsigned int> order_) : type_enum(enum_), order(order_) {};
    };

    explicit Group(GroupType group_type, std::vector<Permutation> generators);

    size_t size_of_permutations() const;

    bool do_groups_commute(const Group& rhs) const;

    std::vector<std::set<size_t>> construct_orbits_of_mults() const;

    bool operator==(const Group& rhs) const;

    bool operator!=(const Group& rhs) const;

    const Group::AlgebraicProperties properties;

    const std::vector<Permutation>& getElements() const;
    const std::vector<Permutation>& getGenerators() const;

    static AlgebraicProperties return_group_info_by_group_type(GroupType group_name);

  private:
    std::vector<Permutation> elements_;
    // There are generators and elements for _specific_ input.
    // Length of Permutation vectors -- number of spins.
    std::vector<Permutation> generators_;
};

struct InitializationError: public std::logic_error {
    explicit InitializationError(const std::string& arg) : logic_error(arg) {}
};

}  // namespace group

group::Group::AlgebraicProperties constructDihedral(unsigned int order);

extern const group::Group::AlgebraicProperties GroupInfoS2;

#endif  //SPINNER_GROUP_H
