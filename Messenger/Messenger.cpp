#include "Messenger.h"
#include <iostream>
#include <cstdlib> // For std::exit

// Static member initialization
std::map<std::string, Messenger> Messenger::instances_;
LogLevel Messenger::globalDefaultLevel_ = LogLevel::WARNING; // change mode
std::mutex Messenger::messengerMutex_;

namespace LogColors {
    const std::string RESET         = "\033[0m";
    const std::string BLACK         = "\033[30m";
    const std::string RED           = "\033[31m";
    const std::string GREEN         = "\033[32m";
    const std::string YELLOW        = "\033[33m";
    const std::string BLUE          = "\033[34m";
    const std::string MAGENTA       = "\033[35m";
    const std::string CYAN          = "\033[36m";
    const std::string WHITE         = "\033[37m";
    // bold
    const std::string BOLD          = "\033[1m";
    const std::string BOLD_RED      = "\033[1;31m";
    const std::string BOLD_YELLOW   = "\033[1;33m";
}

// --- Messenger Implementaion ---
MessengerStream::MessengerStream(Messenger* messenger, LogLevel level, const std::string& source, bool isOnce, int nLimit, const std::string& messageKey)
    : messenger_(messenger), level_(level), source_(source), isOnce_(isOnce), nLimit_(nLimit), messageKey_(messageKey), shouldLog_(true) {
    if(messenger_->getReportLevel() < level_) {
        shouldLog_ = false;
        return;
    }

    if(isOnce_) {
        if(!messenger_->shouldLogOnce(messageKey_)) {
            shouldLog_ = false;
        }
    } else if (nLimit_ > 0) {
        if(!messenger_->shouldLogN(messageKey_, nLimit_)) {
            shouldLog_ = false;
        }
    }
}

MessengerStream::~MessengerStream() {
    if(shouldLog_ && !buffer_.str().empty()) {
        messenger_->print(level_, source_, buffer_.str());
        if(level_ == LogLevel::FATAL) {
            std::exit(EXIT_FAILURE);
        }
    }
}

MessengerStream& MessengerStream::operator<<(std::ostream& (*manip)(std::ostream&)) {
    if(shouldLog_) {
        buffer_ << manip;
    }
    return *this;
}

MessengerStream& MessengerStream::source(const std::string& source) {
    this->source_ = source;
    return *this;
}

// --- Messenger Implementaion ---
Messenger::Messenger() : reportLevel_(globalDefaultLevel_) {}

Messenger::Messenger(LogLevel reportLevel) : reportLevel_(reportLevel) {}

Messenger& Messenger::getInstance(const std::string& instanceName) {
    std::lock_guard<std::mutex> lock(messengerMutex_);
    auto it = instances_.find(instanceName);
    if(it == instances_.end()) {
        it = instances_.emplace(std::piecewise_construct,
                                std::forward_as_tuple(instanceName),
                                std::forward_as_tuple(globalDefaultLevel_)).first;
    }
    return it->second;
}

void Messenger::setDefaultLevel(LogLevel level) {
    std::lock_guard<std::mutex> lock(messengerMutex_);
    globalDefaultLevel_ = level;

    for(auto& pair : instances_) {
        pair.second.setReportLevel(level);
    }
}

void Messenger::setReportLevel(LogLevel level) {
    reportLevel_ = level;
}

LogLevel Messenger::getReportLevel() const {
    return reportLevel_;
}

void Messenger::setSource(const std::string& source) {
    currentSource_ = source;
}

MessengerStream Messenger::fatal(const std::string& source) {
    return MessengerStream(this, LogLevel::FATAL, source.empty() ? currentSource_ : source);
}
MessengerStream Messenger::error(const std::string& source) {
    return MessengerStream(this, LogLevel::ERROR, source.empty() ? currentSource_ : source);
}
MessengerStream Messenger::warning(const std::string& source) {
    return MessengerStream(this, LogLevel::WARNING, source.empty() ? currentSource_ : source);
}
MessengerStream Messenger::status(const std::string& source) {
    return MessengerStream(this, LogLevel::STATUS, source.empty() ? currentSource_ : source);
}
MessengerStream Messenger::info(const std::string& source) {
    return MessengerStream(this, LogLevel::INFO, source.empty() ? currentSource_ : source);
}
MessengerStream Messenger::debug(const std::string& source) {
    return MessengerStream(this, LogLevel::DEBUG, source.empty() ? currentSource_ : source);
}
MessengerStream Messenger::trace(const std::string& source) {
    return MessengerStream(this, LogLevel::TRACE, source.empty() ? currentSource_ : source);
}

// LOG_ONCE variants
MessengerStream Messenger::fatal_once(const std::string& key, const std::string& source) {
    return MessengerStream(this, LogLevel::FATAL, source.empty() ? currentSource_ : source, true, -1, key);
}
MessengerStream Messenger::error_once(const std::string& key, const std::string& source) {
    return MessengerStream(this, LogLevel::ERROR, source.empty() ? currentSource_ : source, true, -1, key);
}
MessengerStream Messenger::warning_once(const std::string& key, const std::string& source) {
    return MessengerStream(this, LogLevel::WARNING, source.empty() ? currentSource_ : source, true, -1, key);
}
MessengerStream Messenger::status_once(const std::string& key, const std::string& source) {
    return MessengerStream(this, LogLevel::STATUS, source.empty() ? currentSource_ : source, true, -1, key);
}
MessengerStream Messenger::info_once(const std::string& key, const std::string& source) {
    return MessengerStream(this, LogLevel::INFO, source.empty() ? currentSource_ : source, true, -1, key);
}


// LOG_N variants
MessengerStream Messenger::fatal_n(const std::string& key, int n, const std::string& source) {
    return MessengerStream(this, LogLevel::FATAL, source.empty() ? currentSource_ : source, false, n, key);
}
MessengerStream Messenger::error_n(const std::string& key, int n, const std::string& source) {
    return MessengerStream(this, LogLevel::ERROR, source.empty() ? currentSource_ : source, false, n, key);
}
MessengerStream Messenger::warning_n(const std::string& key, int n, const std::string& source) {
    return MessengerStream(this, LogLevel::WARNING, source.empty() ? currentSource_ : source, false, n, key);
}
MessengerStream Messenger::status_n(const std::string& key, int n, const std::string& source) {
    return MessengerStream(this, LogLevel::STATUS, source.empty() ? currentSource_ : source, false, n, key);
}
MessengerStream Messenger::info_n(const std::string& key, int n, const std::string& source) {
    return MessengerStream(this, LogLevel::INFO, source.empty() ? currentSource_ : source, false, n, key);
}

bool Messenger::shouldLogOnce(const std::string& key) {
    std::lock_guard<std::mutex> lock(logCountsMutex_);
    if(logCounts_.find(key) == logCounts_.end()) {
        logCounts_[key] = 1; // Mark as logged
        return true;
    }
    return false; // Already logged
}

bool Messenger::shouldLogN(const std::string& key, int n) {
    std::lock_guard<std::mutex> lock(logCountsMutex_);
    if(logCounts_[key] < n) {
        logCounts_[key]++;
        if(logCounts_[key] == n) {
            std::ostringstream suppressMsg;
            suppressMsg << "Message with key [" << key << "] will be suppressed hernceforth (limit " << n << " reached).";
            print(LogLevel::DEBUG, "Messenger", suppressMsg.str());
        }
        return true;
    }
    return false; // Limit reached
}

std::string Messenger::switch_Color(LogLevel level) {
    std::string colorCode;
    switch(level) {
        case LogLevel::FATAL:
            colorCode = LogColors::BOLD_RED;
            break;
        case LogLevel::ERROR:
            colorCode = LogColors::RED;
            break;
        case LogLevel::WARNING:
            colorCode = LogColors::MAGENTA;
            break;
        case LogLevel::STATUS:
            colorCode = LogColors::CYAN;
            break;
        case LogLevel::INFO:
            colorCode = LogColors::GREEN;
            break;
        case LogLevel::DEBUG:
            colorCode = LogColors::YELLOW;
            break;
        case LogLevel::TRACE:
            colorCode = LogColors::BOLD_YELLOW;
            break;
        default:
            colorCode = "";
            break;
    }
    return colorCode;
}

void Messenger::print(LogLevel level, const std::string& source, const std::string& message) {
    if(reportLevel_ < level) {
        return;
    }

    std::ostream& out_stream = (level == LogLevel::ERROR || level == LogLevel::FATAL) ? std::cerr : std::cout;
    std::string levelStr = levelToString(level);
    
    std::string colorCode = switch_Color(level);

    // Simple multi-line handling: indent subsequent lines
    std::istringstream messageStream(message);
    std::string line;
    bool firstLine = true;
    std::string completePrefix;

    while(std::getline(messageStream, line)) {
        if(!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if(firstLine) {
            std::ostringstream prefix_ss;
            prefix_ss << getTimestamp() << colorCode << " [" << levelStr << "] " << LogColors::RESET;
            if(!source.empty()) {
                prefix_ss << "[" << source << "] ";
            }
            completePrefix = prefix_ss.str();
            out_stream << completePrefix << line << std::endl;
            firstLine = false;
        } else {
            if(!completePrefix.empty()) {
                out_stream << std::string(completePrefix.length(), ' ') << line << std::endl;
            } else {
                out_stream << line << std::endl;
            }
        }
    }
}

std::string Messenger::getTimestamp() const {
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    //auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::ostringstream oss;
    oss << LogColors::BOLD << std::put_time(std::localtime(&now_c), "[%Y-%m-%d %H:%M:%S]") << LogColors::RESET;
    // oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
}

std::string Messenger::levelToString(LogLevel level) const {
    switch (level) {
        case LogLevel::FATAL:   return "FATAL"; // Padded for alignment
        case LogLevel::ERROR:   return "ERROR";
        case LogLevel::WARNING: return "WARNING";
        case LogLevel::STATUS:  return "STATUS";
        case LogLevel::INFO:    return "INFO";
        case LogLevel::DEBUG:   return "DEBUG";
        case LogLevel::TRACE:   return "TRACE";
        default:                return "UNKNOWN";
    }
}