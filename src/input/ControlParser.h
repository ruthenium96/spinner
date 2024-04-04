#ifndef SPINNER_CONTROLPARSER_H
#define SPINNER_CONTROLPARSER_H

#include <optional>
#include <yaml-cpp/yaml.h>

#include "src/common/Logger.h"

namespace input {

class ControlParser {
  public:
    explicit ControlParser(YAML::Node control_node);
    common::PrintLevel getPrintLevel() const;
  private:
    std::optional<common::PrintLevel> print_level_;
};

}  // namespace input

#endif  //SPINNER_CONTROLPARSER_H
