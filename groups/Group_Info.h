#ifndef JULY_GROUP_INFO_H
#define JULY_GROUP_INFO_H

#include <cstdint>
#include <vector>

namespace group {
    /*
     group_size : the number of permutations in the group

     number_of_representations : also the number of conjugacy classes

     number_of_generators : the minimal number of elements required to produce all elements of the group

     group_in_form_of_generators : all elements of the group in the form a^vector[0] * b^vector[1] * ...,
                                   where {a, b, ...} -- generators in the _descending_ order of element order.
                                   element order is the minimum m, such that g^m = e
                                   It has the size (group_size) x (number of generators).

     coefficients_of_projectors : has the size (number_of_representations) x (dimension of representation) x (group_size).
                                  important: coefficients_of_projectors[0][0] should correspond to full-symmetric representation.

     */

    struct GroupInfo {
        uint8_t group_size;
        uint8_t number_of_representations;
        std::vector<uint8_t> dimension_of_representation;
        std::vector<uint8_t> number_of_projectors_of_representation;
        uint8_t number_of_generators;
        std::vector<size_t> orders_of_generators;
        std::vector<std::vector<size_t>> group_in_form_of_generators;
        std::vector<size_t> orders_of_elements;
        std::vector<std::vector<std::vector<double>>> coefficients_of_projectors;
    };

    extern const GroupInfo GroupInfoS2;
    extern const GroupInfo GroupInfoS3;

    enum GroupNames {
        S2,
        S3,
    };

    const GroupInfo& return_group_info_by_group_name(GroupNames group_name);
}
#endif //JULY_GROUP_INFO_H
