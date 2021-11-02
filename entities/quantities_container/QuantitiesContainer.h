#ifndef JULY_QUANTITIESCONTAINER_H
#define JULY_QUANTITIESCONTAINER_H

#include "Quantity.h"

#include <map>
#include <vector>

template<class T>
class QuantitiesContainer {
public:
    T& get_unique_quantity_object(const QuantityEnum quantity_enum) {
        if (!contains_quantity(quantity_enum)) {
            throw std::invalid_argument("Quantity " + get_quantity_name(quantity_enum) + " does not exist.");
        }
        size_t index = quantity_enum_to_index_[quantity_enum];
        return unique_quantities_[index];
    }

    const T& get_unique_quantity_object(const QuantityEnum quantity_enum) const {
        if (!contains_quantity(quantity_enum)) {
            throw std::invalid_argument("Quantity " + get_quantity_name(quantity_enum) + " does not exist.");
        }
        size_t index = quantity_enum_to_index_.at(quantity_enum);
        return unique_quantities_[index];
    }

    void emplace_back(T&& object, QuantityEnum quantity_enum) {
        if (contains_quantity(quantity_enum)) {
            throw std::invalid_argument("Quantity " + get_quantity_name(quantity_enum) + " already exists.");
        }
        unique_quantities_.emplace_back(std::move(object));
        quantity_enum_to_index_[quantity_enum] = unique_quantities_.size() - 1;
    }

    void emplace_back(QuantityEnum quantity_enum) {
        if (contains_quantity(quantity_enum)) {
            throw std::invalid_argument("Quantity " + get_quantity_name(quantity_enum) + " already exists.");
        }
        unique_quantities_.emplace_back();
        quantity_enum_to_index_[quantity_enum] = unique_quantities_.size() - 1;
    }

    [[nodiscard]] bool contains_quantity(const QuantityEnum quantity_enum) const {
        return quantity_enum_to_index_.find(quantity_enum) != quantity_enum_to_index_.end();
    }

    using iterator = typename std::map<QuantityEnum, size_t>::iterator;
    iterator begin() { return quantity_enum_to_index_.begin(); }
    iterator end() { return (quantity_enum_to_index_.end()); }
    using const_iterator = typename std::map<QuantityEnum, size_t>::const_iterator;
    [[nodiscard]] const_iterator cbegin() const { return quantity_enum_to_index_.cbegin(); }
    [[nodiscard]] const_iterator cend() const { return (quantity_enum_to_index_.cend()); }

private:
    std::vector<T> unique_quantities_;
    std::map<QuantityEnum, size_t> quantity_enum_to_index_;
};

#endif //JULY_QUANTITIESCONTAINER_H
