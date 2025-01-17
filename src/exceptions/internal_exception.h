#include <exception>
#include <string>

class InternalException : public std::exception {
    private:
        std::string message;

    public:
        explicit InternalException(const std::string &msg);
        const char *what() const noexcept override;
};