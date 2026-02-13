#include "ControlParser.h"

#include "Tools.h"

#ifdef _Eigen_BUILT
    #include "src/linalg_structures/eigen/EigenFactories.h"
#endif
#ifdef _Arma_BUILT
    #include "src/linalg_structures/arma/ArmaFactories.h"
#endif

template<>
std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory> YAML::Node::as() const {
    auto dense_algebra_package_string = as<std::string>();

    if (dense_algebra_package_string == "default") {
        return spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory::defaultFactory();
    } else if (dense_algebra_package_string == "arma") {
#ifdef _Arma_BUILT
        return std::make_shared<spinner::linalg_structures::ArmaDenseTransformAndDiagonalizeFactory>();
#else
        throw std::invalid_argument("Arma was not found, thus cannot be used");
#endif
    } else if (dense_algebra_package_string == "eigen") {
#ifdef _Eigen_BUILT
        return std::make_shared<spinner::linalg_structures::EigenDenseTransformAndDiagonalizeFactory>();
#else
        throw std::invalid_argument("Eigen was not found, thus cannot be used");
#endif
    } else {
        throw std::invalid_argument("Incorrect control::dense_algebra_package: " +
                                    dense_algebra_package_string);
    }
}

namespace spinner::input {
ControlParser::ControlParser(YAML::Node control_node, bool dry_run) {
    print_level_ = extractValue<common::PrintLevel>(control_node, "print_level");
    if (!dry_run) {
        common::Logger::set_level(print_level_.value());
    } else {
        common::Logger::set_level(common::error);
    }

    constructFactoriesList(control_node);

    throw_if_node_is_not_empty(control_node);
}

void ControlParser::constructFactoriesList(YAML::Node& control_node) {
    auto densePrecision =
        extractValue<linalg_structures::Precision>(control_node, "dense_precision");

    auto denseFactory = extractValue<
        std::shared_ptr<linalg_structures::AbstractDenseTransformAndDiagonalizeFactory>
        >(control_node, "dense_algebra_package");
    denseFactory->setPrecision(densePrecision);

    auto sparseFactory =
        linalg_structures::AbstractSparseTransformFactory::defaultSparseFactory();

    factoriesList_ = linalg_structures::FactoriesList(denseFactory, sparseFactory);
}

common::PrintLevel ControlParser::getPrintLevel() const {
    return print_level_.value();
}

const std::optional<linalg_structures::FactoriesList>& ControlParser::getFactoriesList() const {
    return factoriesList_;
}
} // namespace spinner::input