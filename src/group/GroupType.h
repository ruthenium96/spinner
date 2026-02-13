#ifndef SPINNER_GROUPTYPE_H
#define SPINNER_GROUPTYPE_H

#include <optional>

namespace spinner::group {
/**
 * @brief Predefined types of finite groups supported by the Spinner.
 */
enum class GroupTypeEnum {
    S2, ///< Symmetric group of degree 2.
    Dihedral, ///< Dihedral group family (order is specified separately).
};

/**
 * @brief Describes a group by its type and, if necessary, its order.
 */
struct GroupType {
    GroupTypeEnum type_enum; ///< The family or base type of the group.
    std::optional<unsigned int> order; ///< The group order; present if required to specify.

    /**
     * @brief Construct a GroupType from an enumerator only (order left unspecified).
     * @param enum_ The group type enumerator.
     */
    GroupType(GroupTypeEnum enum_) : type_enum(enum_), order(std::nullopt) {};
    /**
     * @brief Construct a GroupType from an enumerator and an optional order.
     * @param enum_   The group type enumerator.
     * @param order_  The group order (may be std::nullopt if order is fixed by the enumerator).
     */
    GroupType(GroupTypeEnum enum_, std::optional<unsigned int> order_) : type_enum(enum_), order(order_) {};
};
} // namespace spinner::group
#endif //SPINNER_GROUPTYPE_H