#ifndef SPINNER_ORDEROFSUMMATION_H
#define SPINNER_ORDEROFSUMMATION_H

#include <memory>
#include <optional>
#include <set>
#include <vector>


namespace index_converter::s_squared {

class OrderOfSummation {
  public:
    struct AdditionInstruction {
        std::vector<uint64_t> positions_of_summands;
        uint64_t position_of_sum;
        // TODO: delete it
        std::optional<size_t> number_of_group = std::nullopt;
    };

    static std::shared_ptr<const OrderOfSummation> constructFromOrbits(
        const std::vector<std::vector<std::set<size_t>>>& all_groups_orbits_of_mults,
        size_t number_of_mults,
        size_t number_of_summation);
    std::vector<AdditionInstruction>::const_iterator begin() const;
    std::vector<AdditionInstruction>::const_iterator end() const;
    const AdditionInstruction& at(size_t i) const;
    const std::vector<AdditionInstruction>& getInstructions() const;
    size_t size() const;

  private:
    std::vector<AdditionInstruction> instructions_;
};
}  // namespace index_converter::s_squared
#endif  //SPINNER_ORDEROFSUMMATION_H
