#include "parameter_exception.h"
#include <string>

ParameterException::ParameterException(const std::string &msg) : message(msg) {}

const char *ParameterException::what() const noexcept {
    return message.c_str();
}
