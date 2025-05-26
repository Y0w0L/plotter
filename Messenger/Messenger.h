#ifndef MESSENGER_H
#define MESSENGER_H

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <map>
#include <mutex> // For thread safety of LOG_ONCE and LOG_N features

class Messenger;

// LOG level
enum class LogLevel {
    FATAL   = 0, // break program
    ERROR   = 1,
    WARNING = 2,
    STATUS  = 3,
    INFO    = 4,
    DEBUG   = 5,
    TRACE   = 6,
    NONE    = 7  // No log
};

// Messenger template class
class MessengerStream {
public:
    MessengerStream(Messenger* messenger,
                    LogLevel level,
                    const std::string& source = "",
                    bool isOnce = false,
                    int nLimit = -1,
                    const std::string& messageKey = "");
    ~MessengerStream();

    template <typename T>
    MessengerStream& operator<<(const T& msg) {
        if(shouldLog_) {
            buffer_ << msg;
        }
        return *this;
    }

    MessengerStream& operator<<(std::ostream& (*manip)(std::ostream&));
    MessengerStream& source(const std::string& source);

private:
    Messenger* messenger_;
    LogLevel level_;
    std::string source_;
    std::ostringstream buffer_;
    bool shouldLog_;

    // For LOG_ONCE and LOG_N
    bool isOnce_;
    int nLimit_;
    std::string messageKey_;
};

class Messenger {
public:
    Messenger();
    explicit Messenger(LogLevel reportLevel);

    static Messenger& getInstance(const std::string& instanceName = "default");
    static void setDefaultLevel(LogLevel level);
    static LogLevel getDefaultLevel();
    static std::string switch_Color(LogLevel level);

    void setReportLevel(LogLevel level);
    LogLevel getReportLevel() const;

    void setSource(const std::string& source);

    MessengerStream fatal(const std::string& source = "");
    MessengerStream error(const std::string& source = "");
    MessengerStream warning(const std::string& source = "");
    MessengerStream status(const std::string& source = "");
    MessengerStream info(const std::string& source = "");
    MessengerStream debug(const std::string& source = "");
    MessengerStream trace(const std::string& source = "");

    // LOG_ONCE and LOG_N variants
    MessengerStream fatal_once(const std::string& key, const std::string& source = "");
    MessengerStream error_once(const std::string& key, const std::string& source = "");
    MessengerStream warning_once(const std::string& key, const std::string& source = "");
    MessengerStream status_once(const std::string& key, const std::string& source = "");
    MessengerStream info_once(const std::string& key, const std::string& source = "");

    MessengerStream fatal_n(const std::string& key, int n, const std::string& source = "");
    MessengerStream error_n(const std::string& key, int n, const std::string& source = "");
    MessengerStream warning_n(const std::string& key, int n, const std::string& source = "");
    MessengerStream status_n(const std::string& key, int n, const std::string& source = "");
    MessengerStream info_n(const std::string& key, int n, const std::string& source = "");

    void print(LogLevel level, const std::string& source, const std::string& message);

    // For LOG_ONCE and LOG_N internal tracking
    bool shouldLogOnce(const std::string& key);
    bool shouldLogN(const std::string& key, int n);

private:
    LogLevel reportLevel_;
    std::string currentSource_;

    static std::map<std::string, Messenger> instances_;
    static LogLevel globalDefaultLevel_;
    static std::mutex messengerMutex_;

    // For LOG_ONCE and LOG_N
    std::map<std::string, int> logCounts_;
    std::mutex logCountsMutex_;

    std::string getTimestamp() const;
    std::string levelToString(LogLevel level) const;
};

#define LOG_KEY (std::string(__FILE__) + ":" + std::to_string(__LINE__))

#define MSG_FATAL Messenger::getInstance().fatal_once(LOG_KEY)
#define MSG_ERROR Messenger::getInstance().error_once(LOG_KEY)
#define MSG_WARNING Messenger::getInstance().warning_once(LOG_KEY)
#define MSG_STATUS Messenger::getInstance().status_once(LOG_KEY)
#define MSG_INFO Messenger::getInstance().info_once(LOG_KEY)

#define LOG_FATAL Messenger::getInstance().fatal()
#define LOG_ERROR Messenger::getInstance().error()
#define LOG_WARNING Messenger::getInstance().warning()
#define LOG_STATUS Messenger::getInstance().status()
#define LOG_INFO Messenger::getInstance().info()
#define LOG_DEBUG Messenger::getInstance().debug()
#define LOG_TRACE Messenger::getInstance().trace()

// LOG_N variants
#define LOG_FATAL_N(N) Messenger::getInstance().fatal_n(LOG_KEY, N)
#define LOG_ERROR_N(N) Messenger::getInstance().error_n(LOG_KEY, N)
#define LOG_WARNING_N(N) Messenger::getInstance().warning_n(LOG_KEY, N)
#define LOG_STATUS_N(N) Messenger::getInstance().status_n(LOG_KEY, N)
#define LOG_INFO_N(N) Messenger::getInstance().info_n(LOG_KEY, N)

#endif // MESSENGER_H