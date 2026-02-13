#ifndef SPINNER_CONTROLPARSER_H
#define SPINNER_CONTROLPARSER_H

#include <optional>
#include <yaml-cpp/yaml.h>

#include "src/common/Logger.h"
#include "src/linalg_structures/FactoriesList.h"

namespace spinner::input {

class ControlParser {
  public:
    ControlParser(YAML::Node control_node, bool dry_run);
    common::PrintLevel getPrintLevel() const;
    const std::optional<linalg_structures::FactoriesList>& getFactoriesList() const;
  private:
    void constructFactoriesList(YAML::Node& control_node);

    std::optional<common::PrintLevel> print_level_;
    std::optional<linalg_structures::FactoriesList> factoriesList_;
};

} // namespace spinner::input

#endif  //SPINNER_CONTROLPARSER_H
