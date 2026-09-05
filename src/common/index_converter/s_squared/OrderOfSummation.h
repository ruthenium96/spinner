#ifndef SPINNER_ORDEROFSUMMATION_H
#define SPINNER_ORDEROFSUMMATION_H

#include <cstddef>
#include <memory>
#include <set>
#include <vector>


namespace spinner::index_converter::s_squared {

class OrderOfSummation {
  public:
    struct AdditionInstruction {
        std::size_t left;
        std::size_t right;
        std::size_t result;
    };

    static std::shared_ptr<const OrderOfSummation> constructFromOrbits(
        const std::vector<std::vector<std::set<size_t>>>& all_groups_orbits_of_mults,
        size_t number_of_mults);
    std::vector<AdditionInstruction>::const_iterator begin() const;
    std::vector<AdditionInstruction>::const_iterator end() const;
    const AdditionInstruction& at(size_t i) const;
    const std::vector<AdditionInstruction>& getInstructions() const;
    size_t size() const;

  private:
    std::vector<AdditionInstruction> instructions_;
};
} // namespace spinner::index_converter::s_squared
#endif  //SPINNER_ORDEROFSUMMATION_H
