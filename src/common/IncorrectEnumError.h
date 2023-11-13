#ifndef SPINNER_INCORRECTENUMERROR_H
#define SPINNER_INCORRECTENUMERROR_H

#include <stdexcept>

namespace common {

class IncorrectEnumError: public std::logic_error {
  public:
    explicit IncorrectEnumError(const std::string& arg) : logic_error(arg) {}
};

}  // namespace common

#endif  //SPINNER_INCORRECTENUMERROR_H
