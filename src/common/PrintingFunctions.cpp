#include "PrintingFunctions.h"

std::ostream& operator<<(std::ostream& os, const space::Space& space) {
    for (const space::Subspace& subspace : space.getBlocks()) {
        os << subspace;
    }
    os << "------" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const space::Subspace& subspace) {
    os << subspace.properties;
    subspace.decomposition->print(os);
    if (subspace.dense_semiunitary_matrix.has_value()) {
        os << std::endl;
        subspace.dense_semiunitary_matrix.value()->print(os);
    }
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Spectrum& spectrum) {
    for (const Subspectrum& subspectrum : spectrum.blocks) {
        os << subspectrum;
    }
    os << "------" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum) {
    os << subspectrum.properties << std::endl;
    subspectrum.raw_data->print(os);
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (const Submatrix& submatrix : matrix.blocks) {
        os << submatrix;
    }
    os << "------" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix) {
    os << submatrix.properties;
    submatrix.raw_data->print(os);
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const BlockProperties& properties) {
    os << "Total n-projection: ";
    if (properties.n_proj.has_value()) {
        os << properties.n_proj.value();
    } else {
        os << "none";
    }
    os << "\nTotal multiplicity: ";
    if (properties.total_mult.has_value()) {
        os << properties.total_mult.value();
    } else {
        os << "none";
    }
    os << '\n'
       << "    dimensionality: " << properties.dimensionality << '\n'
       << "        degeneracy: " << properties.degeneracy << '\n'
       << "    representation: " << properties.get_representation_name() << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const common::QuantityEnum& quantity_enum) {
    os << "Quantity type: " << common::get_quantity_name(quantity_enum) << std::endl;
    return os;
}