#pragma once

namespace entities {
class Entity {
  public:
    enum Type { SPACE, OPERATOR, MATRIX, SPECTRUM };
    struct History {
        bool isTzSorted = false;
        bool isC2Symmetrized = false;
    };
    Entity(Type t) : type(t) {
    }
    Type type;
    History history;
};
} // namespace entities