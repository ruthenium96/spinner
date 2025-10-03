#include "tests/integration_tests/spectrum_equivalence/equivalence_check.h"

#include <algorithm>
#include <cstddef>
#include <optional>
#include "magic_enum.hpp"

#include "src/common/OneOrMany.h"
#include "src/common/Quantity.h"

using QuantumValues = std::array<std::optional<double>, magic_enum::enum_count<common::QuantityEnum>()>;

inline std::vector<QuantumValues> construct_final_vector(runner::Runner& runner) {
    std::vector<QuantumValues> vector;

    auto factory = runner.getDataStructuresFactories();
    auto degeneracy_vector = factory.createVector();

    std::array<std::optional<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>, 
        magic_enum::enum_count<common::QuantityEnum>()> values_vectors;

    for (const auto& quantity_enum_ : magic_enum::enum_values<common::QuantityEnum>()) {
        if (runner.getSpectrum(quantity_enum_).has_value()) {
            values_vectors[magic_enum::enum_integer<common::QuantityEnum>(quantity_enum_)] 
                = factory.createVector();
        } else {
            values_vectors[magic_enum::enum_integer<common::QuantityEnum>(quantity_enum_)] 
                = std::nullopt;
        }
    }

    EXPECT_TRUE(holdsOne(runner.getSpectrum(common::Energy).value()));
    SpectrumRef energy_spectrum_ref = getOneRef(runner.getSpectrum(common::Energy).value());
    for (const auto& subspectrum_ref : energy_spectrum_ref.blocks) {
        const auto& subspectrum = subspectrum_ref.get();
        degeneracy_vector->add_identical_values(
            subspectrum.raw_data->size(),
            subspectrum.properties.degeneracy);
    }

    for (int j = 0; j < values_vectors.size(); ++j) {
        const auto& quantity_enum_ = magic_enum::enum_value<common::QuantityEnum>(j);
        if (values_vectors[j].has_value()) {
            SpectrumRef quantity_spectrum_ref = getOneRef(runner.getSpectrum(quantity_enum_).value());
            for (const auto& subspectrum_ref : quantity_spectrum_ref.blocks) {
                const auto& subspectrum = subspectrum_ref.get();
                values_vectors[j].value().get()->concatenate_with(subspectrum.raw_data);
            }
        }
    }

    values_vectors[magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value().get()->subtract_minimum();

    for (size_t i = 0; i < degeneracy_vector->size(); ++i) {
        QuantumValues quantum_values;
        for (size_t j = 0; j < values_vectors.size(); ++j) {
            if (values_vectors[j].has_value()) {
                quantum_values[j] = values_vectors[j].value()->at(i);
            }
        }
        for (size_t j = 0; j < (size_t)degeneracy_vector->at(i); ++j) {
            vector.push_back(quantum_values);
        }
    }
    std::stable_sort(vector.begin(), vector.end());
    return vector;
}

void expect_final_vectors_equivalence(runner::Runner& simple, runner::Runner& second, bool assert) {

    auto first_vector = construct_final_vector(simple);
    auto second_vector = construct_final_vector(second);

    ASSERT_EQ(first_vector.size(), second_vector.size());

    QuantumValues quantum_values_sum_first;
    QuantumValues quantum_values_sum_second;

    for (int j = 0; j < quantum_values_sum_first.size(); ++j) {
        if (j == magic_enum::enum_integer<common::QuantityEnum>(common::Energy)) {
            continue;
        }
        if (first_vector.at(0)[j].has_value()) {
            quantum_values_sum_first[j] = 0;
        }
        if (second_vector.at(0)[j].has_value()) {
            quantum_values_sum_second[j] = 0;
        }
    }

    double last_energy = INFINITY;

    for (size_t i = 0; i < first_vector.size(); ++i) {
        // TODO: epsilon
        if (assert) {
            ASSERT_NEAR(first_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value(), 
                second_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value(), 1e-9);
        } else {
            EXPECT_NEAR(first_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value(), 
                second_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value(), 1e-9);
        }

        // TODO: values may differ, but the sums of values for degenerate eigenvectors should be the same
        if (std::abs(last_energy - first_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value()) > 1e-5) {
            for (int j = 0; j < quantum_values_sum_first.size(); ++j) {
                for (int k = 0; k < quantum_values_sum_second.size(); ++k) {
                    if (quantum_values_sum_first[j].has_value() 
                        && quantum_values_sum_second[k].has_value()) {
                        if (magic_enum::enum_value<common::QuantityEnum>(j) != common::S_total_squared) {
                            quantum_values_sum_first[j].value() *= 3;
                        }
                        if (magic_enum::enum_value<common::QuantityEnum>(k) != common::S_total_squared) {
                            quantum_values_sum_second[k].value() *= 3;
                        }
                        if (assert) {
                            ASSERT_NEAR(quantum_values_sum_first[j].value(), quantum_values_sum_second[k].value(), 1e-3)
                            << "Names of values, j: " << magic_enum::enum_entries<common::QuantityEnum>()[j].second
                            << ", k: " << magic_enum::enum_entries<common::QuantityEnum>()[k].second;
                        } else {
                            EXPECT_NEAR(quantum_values_sum_first[j].value(), quantum_values_sum_second[k].value(), 1e-3)
                            << "Names of values, j: " << magic_enum::enum_entries<common::QuantityEnum>()[j].second
                            << ", k: " << magic_enum::enum_entries<common::QuantityEnum>()[k].second;
                        }
                        quantum_values_sum_first[j] = 0;
                        quantum_values_sum_second[k] = 0;
                    }
                }
            }
        }
        for (int j = 0; j < quantum_values_sum_first.size(); ++j) {
            const auto& quantity_enum_ = magic_enum::enum_value<common::QuantityEnum>(j);
            if (quantum_values_sum_first[j].has_value()) {
                quantum_values_sum_first[j].value() += first_vector[i][quantity_enum_].value();
            }
            if (quantum_values_sum_second[j].has_value()) {
                quantum_values_sum_second[j].value() += second_vector[i][quantity_enum_].value();  
            }
        }
        last_energy = first_vector[i][magic_enum::enum_integer<common::QuantityEnum>(common::Energy)].value();
    }
}

std::string spectrum_final_equivalence_test_name_generator(const ::testing::TestParamInfo<TestParam>& info) {
    const auto& [optimization_list, mults] = info.param;
    std::string name; 

    name += std::string(optimization_list.isLexBasis() ? "lex" : "ito");
    if (optimization_list.isTzSorted()) {
        name += "_tz";
    }
    if (optimization_list.isTSquaredSorted()) {
        name += "_t2";
    }
    if (optimization_list.isPositiveProjectionsEliminated()) {
        name += "_posElim";
    }
    if (optimization_list.isNonMinimalProjectionsEliminated()) {
        name += "_nonminElim";
    }
    if (optimization_list.isSSquaredTransformed()) {
        name += "_s2";
    }
    if (optimization_list.isFTLMApproximated()) {
        name += "_FTLM";
    }
    if (!optimization_list.getGroupsToApply().empty()) {
        for (int i = 0; i < optimization_list.getGroupsToApply().size(); ++i) {
            name += "_g" + std::to_string(i) + "_";
            const auto& group = optimization_list.getGroupsToApply()[i];
            for (int j = 0; j < group.getGenerators().size(); ++j) { // skip first
                const auto& gen = group.getGenerators()[j];                
                for (int k = 0; k < gen.size(); ++k) {
                    name += std::to_string(gen[k]);
                    if (k + 1 < gen.size()) {
                        name += "x";
                    }
                }
                if (j + 1 < group.getGenerators().size()) {
                    name += "_";
                }
            }
        }
        if (optimization_list.isNonAbelianSimplified()) {
            name += "_nonabSimp";
        }
    }
    name += "_";
    for (int i = 0; i < mults.size(); ++i) {
        name += std::to_string(mults[i]);
        if (i + 1 < mults.size()) {
            name += "x";
        }
    }

    return name;
}

std::ostream& common::physical_optimization::operator<<(std::ostream& os, const common::physical_optimization::OptimizationList& optimization_list) {
    os << std::string(optimization_list.isLexBasis() ? "LEX|" : "ITO|"); 
    if (optimization_list.isTzSorted()) { 
        os << "tz|"; 
    } 
    if (optimization_list.isTSquaredSorted()) { 
        os << "t^2|"; 
    } 
    if (optimization_list.isPositiveProjectionsEliminated()) { 
        os << "posElim|"; 
    } 
    if (optimization_list.isNonMinimalProjectionsEliminated()) { 
        os << "nonminElim|"; 
    } 
    if (optimization_list.isSSquaredTransformed()) { 
        os << "s^2|"; 
    } 
    if (optimization_list.isFTLMApproximated()) { 
        os << "FTLM|"; 
    } 
    if (!optimization_list.getGroupsToApply().empty()) {
        os << "Groups:{"; 
        for (int i = 0; i < optimization_list.getGroupsToApply().size(); ++i) { 
            os << "["; 
            const auto& group = optimization_list.getGroupsToApply()[i]; 
            for (int j = 0; j < group.getGenerators().size(); ++j) { 
                const auto& el = group.getGenerators()[j];
                os << "(";
                for (int k = 0; k < el.size(); ++k) { 
                    os << std::to_string(el[k]);
                    if (k + 1 < el.size()) {
                        os << ", ";
                    }
                } 
                os << ")"; 
                if (j + 1 < group.getGenerators().size()) {
                    os << ", ";
                }
            }
            os << "]";
            if (i + 1 < optimization_list.getGroupsToApply().size()) {
                os << "; ";
            }
        } 
        os << "}";
        if (optimization_list.isNonAbelianSimplified()) { 
            os << "|nonabSimp"; 
        } 
    }
    return os;
}