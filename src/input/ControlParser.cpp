#include "ControlParser.h"

#include "Tools.h"

namespace input {
ControlParser::ControlParser(YAML::Node control_node) {
    print_level_ = extractValue<common::PrintLevel>(control_node, "print_level");
    common::Logger::set_level(print_level_.value());

    auto densePrecision =
        extractValue<quantum::linear_algebra::Precision>(control_node, "dense_precision");

    auto denseFactory =
        quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory::defaultFactory();
    denseFactory->setPrecision(densePrecision);

    auto sparseFactory =
        quantum::linear_algebra::AbstractSparseTransformFactory::defaultSparseFactory();

    factoriesList_ = quantum::linear_algebra::FactoriesList(denseFactory, sparseFactory);

    throw_if_node_is_not_empty(control_node);
}

common::PrintLevel ControlParser::getPrintLevel() const {
    return print_level_.value();
}

const std::optional<quantum::linear_algebra::FactoriesList>& ControlParser::getFactoriesList() const {
    return factoriesList_;
}
}  // namespace input