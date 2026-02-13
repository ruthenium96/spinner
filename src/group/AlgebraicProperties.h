#ifndef SPINNER_ALGEBRAICPROPERTIES_H
#define SPINNER_ALGEBRAICPROPERTIES_H

#include <cstdint>
#include <vector>

#include "GroupType.h"

namespace spinner::group {

/**
 * @struct AlgebraicProperties
 * @brief Contains key algebraic properties of a finite group @f$G@f$, its matrix representations 
 * and all their projectors.
 */
struct AlgebraicProperties {
    /// @brief The total number of elements in the group.
    uint8_t group_size;
    /// @brief The number of irreducible representations / conjugacy classes of the group.
    uint8_t number_of_representations;
    /// @brief Indicates if the group is Abelian (commutative).
    bool is_abelian;
    /**
     * @brief Dimensions of each irreducible representation. 
     * `dimensions_of_representations[i]` is the dimension of the i-th
     * irreducible representation. Its size is `number_of_representations`.
     */
    std::vector<uint8_t> dimensions_of_representations;
    /**
     * @brief The number of projectors for each representation, 
     * squared dimension of representation. Size is `number_of_representations`.
     */
    std::vector<uint8_t> number_of_projectors_of_representation;
    /// @brief The minimum number of generating elements required to build the group.
    uint8_t number_of_generators;
    /**
     * @brief `orders_of_generators[j]` is the smallest positive integer @f$m@f$ such that
     * @f$g_j^m = e@f$, where @f$g_j@f$ is the j-th generator. Size is `number_of_generators`.
     * @note Generators should be stored in **descending order of their element
     * order** (i.e., `orders_of_generators[0]` is the largest).
     */
    std::vector<size_t> orders_of_generators;
    /**
     * @brief Complete enumeration of group elements expressed
     * as product of powers of generators. A 2D vector of size `group_size` x `number_of_generators`. 
     * An order of multiplication of generators defined in @ref group::Group "the constructor of Group".
     */
    std::vector<std::vector<size_t>> group_in_form_of_generators;
    /**
     * @brief `orders_of_elements[j]` is the smallest positive integer @f$m@f$ such that
     * @f$g_j^m = e@f$, where @f$g_j@f$ is the j-th element of the group. The indexing of elements
     *  matches the order of `group_in_form_of_generators`. Size is `group_size`.
     */
    std::vector<size_t> orders_of_elements;
    /**
     * @brief Coefficients for constructing projection operators onto irreducible subspaces.
     * It has the size:
     * `number_of_representations` x `dimensions_of_representations[i]` x `group_size`.
     * `coefficients_of_projectors[irrep_idx][component][element_idx]` provides the
     * coefficient (derived from the representation matrix) needed to project
     * onto the `component`-th basis vector of the `irrep_idx`-th
     * irreducible representation, using the `element_idx`-th group element in the sum.
     * @note The first representation (`irrep_idx = 0`) **must** correspond
     * to the trivial (full-symmetric) representation. Therefore,
     * `coefficients_of_projectors[0][0]={1, 1, ..., 1}`. 
     */
    std::vector<std::vector<std::vector<double>>> coefficients_of_projectors;
    
    /**
     * @brief Static contructor of the AlgebraicProperties.
     * 
     * @param group_type 
     * @return const AlgebraicProperties 
     */
    static const AlgebraicProperties constructAlgebraicProperties(GroupType group_type);
private:
    AlgebraicProperties() = default;
};

} // namespace spinner::group
#endif  //SPINNER_ALGEBRAICPROPERTIES_H
