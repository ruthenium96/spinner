#ifndef JULY_GROUP_H
#define JULY_GROUP_H

#include <cstdint>
#include <vector>

using Permutation = std::vector<uint8_t>;

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
    };

    enum GroupTypeEnum {
        S2,
        S3,
    };

    explicit Group(GroupTypeEnum group_name, std::vector<Permutation> generators);

    [[nodiscard]] std::vector<std::vector<uint8_t>>
    permutate(const std::vector<uint8_t>& initial) const;

    bool operator==(const Group& rhs) const;

    bool operator!=(const Group& rhs) const;

    const Group::AlgebraicProperties& properties;

    std::vector<Permutation> elements_;

    static const AlgebraicProperties& return_group_info_by_group_name(GroupTypeEnum group_name);

  private:
    // There are generators and elements for _specific_ input.
    // Length of Permutation vectors -- number of spins.
    std::vector<Permutation> generators_;
};

extern const Group::AlgebraicProperties GroupInfoS2;
extern const Group::AlgebraicProperties GroupInfoS3;

#endif  //JULY_GROUP_H
