#ifndef __DCA_LOGGER_H__
#define __DCA_LOGGER_H__

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_DEFAULT "\x1b[0m"

#include <iostream>
#include <sstream>
#include <string>

#define LOG DCA::LogMessage(__FILE__, __func__, __LINE__)

namespace DCA {

class Logger {
public:
    void logIt(const std::string& s);
};

class LogMessage {
public:
    LogMessage(const char* file, const char* function, int line);

    template <typename T>
    inline LogMessage& operator<<(const T& t) {
        os << t;
        return *this;
    }

    ~LogMessage();

private:
    std::ostringstream os;
};

}  // namespace DCA

#endif /* __DCA_LOGGER_H__ */