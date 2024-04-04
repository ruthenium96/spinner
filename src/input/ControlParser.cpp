#include "ControlParser.h"

#include "Tools.h"

namespace input {
ControlParser::ControlParser(YAML::Node control_node) {
    print_level_ = extractValue<common::PrintLevel>(control_node, "print_level");
    common::Logger::set_level(print_level_.value());

    throw_if_node_is_not_empty(control_node);
}

common::PrintLevel ControlParser::getPrintLevel() const {
    return print_level_.value();
}
}  // namespace input