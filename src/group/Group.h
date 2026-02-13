#ifndef SPINNER_GROUP_H
#define SPINNER_GROUP_H

#include <set>
#include <vector>

#include "AlgebraicProperties.h"

namespace spinner::group {

using PermutationArrayType = uint8_t;
using Permutation = std::vector<PermutationArrayType>;

/**
 * @class Group
 * @brief Represents a finite group defined by permutation generators.
 *
 * This class encapsulates a group of permutations, its algebraic properties,
 * and provides operations such as group generation, comparison, orbit computation,
 * and commutativity checks. The group is fully determined by its generator set
 * and the intended group type (which determines its expected structure).
 *
 * Upon construction, the class:
 *   - Verifies that the provided generators meet the invariants (size, validity,
 *     order properties) expected from the algebraic properties.
 *   - Explicitly constructs all group elements by taking products of the generators
 *     according to a canonical word representation.
 *   - Stores the computed algebraic properties (size, conjugacy classes, character data,
 *     projector coefficients, etc.) in a read-only member.
 *
 * Two groups are considered equal if they are isomorphic. 
 */
class Group {
  public:
    /**
     * @brief Constructs a group from its type and a set of permutation generators.
     *
     * @param group_type  The intended type of the group (e.g., S2, Dihedral).
     * @param generators  A non‑empty vector of permutations that generate the group.
     *                    All permutations must have the same size (number of multiplcities)
     *                    and be valid bijections on the set {0,1,...,n-1}.
     *
     * @throws InitializationError if any of the following checks fail:
     *   - The number of generators does not match the expected value from `group_type`.
     *   - Generators have inconsistent or invalid sizes.
     *   - A generator raised to its expected order does not yield the identity.
     *   - A generator raised to a power smaller than its order already equals the identity.
     *   - After constructing all elements, an element’s order is inconsistent.
     *
     * @note The generators are assumed to be stored internally in the order
     *       required by the algebraic properties (descending order of element order).
     */
    explicit Group(GroupType group_type, std::vector<Permutation> generators);

    /**
     * @brief Returns the number of multiplicities on which the permutations act.
     * @return The size (length) of each permutation vector.
     */
    size_t size_of_permutations() const;

    /**
     * @brief Checks whether every generator of this group commutes with every
     *        generator of another group.
     *
     * @param rhs The other group to compare with.
     * @return true  if all pairs of generators commute,
     * @return false otherwise.
     */
    bool do_groups_commute(const Group& rhs) const;

    /**
     * @brief Computes the orbits of the group action on the set of multiplicity centers.
     *
     * @return std::vector<std::set<size_t>>
     *         A vector of disjoint orbits; the order of orbits is unspecified.
     */
    std::vector<std::set<size_t>> construct_orbits_of_mults() const;

    /**
     * @brief Compares two groups for equality.
     *
     * @details Two groups are considered equal if they isomorphic. The actual generators or the
     *          internal order of elements are irrelevant for equality.
     *
     * @param rhs The other group to compare with.
     * @return true  if both groups contain exactly the same set of permutations,
     * @return false otherwise.
     */
    bool operator==(const Group& rhs) const;

    /**
     * @brief Compares two groups for inequality.
     *
     * @param rhs The other group to compare with.
     * @return true  if the groups are not isomorphic,
     * @return false otherwise.
     */
    bool operator!=(const Group& rhs) const;

    /**
     * @brief Algebraic properties of the group (size, representation data, etc.).
     *
     * This member is immutable after construction and is guaranteed to be
     * consistent with the actual permutation elements stored in the group.
     */
    const AlgebraicProperties properties;

    /**
     * @brief Provides read‑only access to the full list of group elements.
     * @return A const reference to the internal vector of permutations,
     *         indexed in the order determined by 
     *         @ref group::AlgebraicProperties "AlgebraicProperties".
     */
    const std::vector<Permutation>& getElements() const;
    /**
     * @brief Provides read‑only access to the generators of the group.
     * @return A const reference to the internal vector of generator permutations.
     */
    const std::vector<Permutation>& getGenerators() const;

  private:
    std::vector<Permutation> elements_;
    std::vector<Permutation> generators_;
};

} // namespace spinner::group

#endif  //SPINNER_GROUP_H
