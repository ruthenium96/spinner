#include "FTLMImplicitQuantityEigendecompositor.h"
#include "src/common/Quantity.h"
#include "src/eigendecompositor/ImplicitQuantityEigendecompositor.h"

namespace eigendecompositor {

FTLMImplicitQuantityEigendecompositor::FTLMImplicitQuantityEigendecompositor(
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
    quantum::linear_algebra::FactoriesList factories_list,
    common::QuantityEnum quantity_implicit_enum,
    uint32_t max_ntz_proj) : 
    ImplicitQuantityEigendecompositor(
        std::move(eigendecompositor),
        factories_list, 
        quantity_implicit_enum,
        max_ntz_proj),
    quantity_implicit_enum_(quantity_implicit_enum) {}

std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
FTLMImplicitQuantityEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {

    auto mb_unitary_transformation_matrix = ImplicitQuantityEigendecompositor::BuildSubspectra(
        number_of_block,
        subspace);

    const auto& quantity_implicit_enum_subspectrum = ImplicitQuantityEigendecompositor::getSubspectrum(quantity_implicit_enum_, number_of_block).value().get();
    const auto& squared_back_projection_subspectrum = getSubspectrum(common::squared_back_projection, number_of_block).value().get();

    auto raw_data = quantity_implicit_enum_subspectrum.raw_data->element_wise_multiplication(squared_back_projection_subspectrum.raw_data);
    auto properties = subspace.properties;
    properties.degeneracy *= (double)subspace.size();

    quantity_implicit_spectrum_multiplied_by_squared_back_projection_.blocks[number_of_block] =
        Subspectrum(std::move(raw_data), properties);

    return mb_unitary_transformation_matrix;
}

std::optional<SpectrumRef>
FTLMImplicitQuantityEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == quantity_implicit_enum_) {
        return SpectrumRef(quantity_implicit_spectrum_multiplied_by_squared_back_projection_);
    }
    return ImplicitQuantityEigendecompositor::getSpectrum(quantity_enum);
}

std::optional<std::reference_wrapper<const Subspectrum>>
FTLMImplicitQuantityEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == quantity_implicit_enum_) {
        return quantity_implicit_spectrum_multiplied_by_squared_back_projection_.blocks[number_of_block];
    }
    return ImplicitQuantityEigendecompositor::getSubspectrum(quantity_enum, number_of_block);
}

void FTLMImplicitQuantityEigendecompositor::initialize(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    uint32_t number_of_subspaces) {
    if (quantity_implicit_spectrum_multiplied_by_squared_back_projection_.blocks.empty()) {
        quantity_implicit_spectrum_multiplied_by_squared_back_projection_.blocks.resize(number_of_subspaces);
    }
    ImplicitQuantityEigendecompositor::initialize(
        operators_to_calculate,
        derivatives_operators_to_calculate,
        number_of_subspaces);
}

void FTLMImplicitQuantityEigendecompositor::finalize() {
    ImplicitQuantityEigendecompositor::finalize();
};

} // namespace eigendecompositor