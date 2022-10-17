#ifndef SPINNER_LOGGER_H
#define SPINNER_LOGGER_H

#include <ctime>
#include <iostream>
#include <mutex>

extern std::mutex logMutex;
// clang-format off
#define LOG(t,x) do{\
    const std::lock_guard<std::mutex> lock(logMutex);\
    time_t now = time(0);\
    char timestamp[10] = "";\
    strftime (timestamp, 10,"%H:%M:%S", localtime(&now));\
    std::cout << "> " << timestamp << " " << t << " " << x << std::endl;\
}while(false);
#define LOG_DEBUG(x)   LOG("DBG:",x);
#define LOG_INFO(x)    LOG("\033[21;34mINF:\033[0m", x);
#define LOG_WARNING(x) LOG("\033[21;33mWRN:\033[0m", x);
#define LOG_ERROR(x)   LOG("\033[1;31mERR:\033[0m", x);
// clang-format on
// for color text
// https://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal

#endif  //SPINNER_LOGGER_H