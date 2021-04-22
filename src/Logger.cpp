#include <DCA/Logger.h>

namespace DCA {

void Logger::logIt(const std::string& s) {
    std::cout << s << ANSI_COLOR_DEFAULT << std::endl;
}

LogMessage::LogMessage(const char* file, const char* function, int line) {
    os << file << ": " << function << "(" << line << ") ";
}

LogMessage::~LogMessage() {
    Logger logger;
    logger.logIt(os.str());
}

}  // namespace DCA