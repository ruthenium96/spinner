#ifndef SPINNER_GROUPTYPE_H
#define SPINNER_GROUPTYPE_H

#include <optional>

namespace spinner::group {
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
} // namespace spinner::group
#endif //SPINNER_GROUPTYPE_H