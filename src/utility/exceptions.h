#pragma once
#include <stdexcept>

class BaseException : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class ParameterException final : public BaseException {
public:
    using BaseException::BaseException;
};

class FileException final : public BaseException {
public:
    using BaseException::BaseException;
};

class InternalException final : public BaseException {
public:
    using BaseException::BaseException;
};
