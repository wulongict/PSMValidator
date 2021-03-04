//
// Created by wulong on 2/1/18.
//

#include "CDebugMode.h"
#include "spdlog/sinks/sink.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"
#include "spdlog/sinks/ansicolor_sink.h"

using namespace std;
using namespace spdlog;

CDebugMode * CDebugMode::pDebug = nullptr;

CDebugMode *CDebugMode::callDebug() {
    if(CDebugMode::pDebug == nullptr)
    {
        CDebugMode::pDebug = new CDebugMode();

    }
    return CDebugMode::pDebug;

}

void CDebugMode::releasePtr() {
    if(pDebug!=nullptr) {
        delete pDebug;
        pDebug = nullptr;
    }
}


void initlog(string logfile, string logname) {
    auto filelog = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile.c_str(),
                                                                          1024*1024*10, 3);
    auto stdsink = std::make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>();
    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(filelog);
    sinks.push_back(stdsink);
    auto combined_logger = std::make_shared<spdlog::logger>(logname, std::begin(sinks), std::end(sinks));
    spdlog::register_logger(combined_logger);

}
