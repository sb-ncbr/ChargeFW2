#include <string>
#include "internal_exception.h"

InternalException::InternalException(const std::string &msg) : message(msg) {}

const char *InternalException::what() const noexcept {
    return message.c_str();
}
