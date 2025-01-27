#include <string>
#include "file_exception.h"

FileException::FileException(const std::string &msg) : message(msg) {}

const char *FileException::what() const noexcept {
    return message.c_str();
}
