#pragma once

#include <exception>
#include <string>

class ParameterException : public std::exception {
    private:
        std::string message;

    public:
        explicit ParameterException(const std::string &msg);
        const char *what() const noexcept override;
};
