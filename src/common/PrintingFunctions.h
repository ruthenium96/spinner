#ifndef SPINNER_PRINTINGFUNCTIONS_H
#define SPINNER_PRINTINGFUNCTIONS_H

#include <ostream>

#include "Quantity.h"
#include "src/entities/BlockProperties.h"
#include "src/entities/matrix/Matrix.h"
#include "src/entities/matrix/Submatrix.h"
#include "src/entities/spectrum/Spectrum.h"
#include "src/entities/spectrum/Subspectrum.h"
#include "src/space/Space.h"
#include "src/space/Subspace.h"

std::ostream& operator<<(std::ostream& os, const space::Space& space);
std::ostream& operator<<(std::ostream& os, const space::Subspace& subspace);
std::ostream& operator<<(std::ostream& os, const Spectrum& spectrum);
std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum);
std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix);
std::ostream& operator<<(std::ostream& os, const BlockProperties& properties);
std::ostream& operator<<(std::ostream& os, const common::QuantityEnum& quantity);


#endif  //SPINNER_PRINTINGFUNCTIONS_H
