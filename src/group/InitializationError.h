#ifndef SPINNER_GROUP_INITIALIZATIONERROR_H
#define SPINNER_GROUP_INITIALIZATIONERROR_H

#include <stdexcept>
#include <string>

namespace spinner::group {
/**
 * @brief Exception thrown when group initialization fails.
 */
struct InitializationError: public std::logic_error {
    explicit InitializationError(const std::string& arg) : logic_error(arg) {}
};
} // namespace spinner::group
#endif //SPINNER_GROUP_INITIALIZATIONERROR_H