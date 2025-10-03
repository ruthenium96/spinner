#include "src/common/Logger.h"
#include "src/common/runner/Executer.h"
#include "src/input/Parser.h"

const std::string help_string = "\nUsage: ./spinner_main filename.yml [--dry]";

int main(int argc, char *argv[]) {
    common::Logger::set_pattern("[%C/%m/%d %H:%M:%S.%e] [%t] %v");
    if (argc == 1) {
        common::Logger::error_msg("Missing input file.{}", help_string);
        return -1;
    }
    if (argc > 3) {
        common::Logger::error_msg("Too much arguments.{}", help_string);
        return -1;
    }
    bool dry_run = false;
    if (argc == 3 && argv[2] != std::string("--dry")) {
        common::Logger::error_msg("Unknown argument: {}{}", argv[2], help_string);
        return -1;
    }
    if (argc == 3 && argv[2] == std::string("--dry")) {
        dry_run = true;
    }
    try {
        if (dry_run) {
            auto parser = input::Parser(argv[1], dry_run);
            runner::Executer::dry_execute(parser);
        } else {
            common::Logger::basic_msg("The calculations have started");
            common::Logger::separate(0, common::basic);        
            auto parser = input::Parser(argv[1]);
            runner::Executer::execute(parser);
        }
    } catch (const std::exception& e) {
        common::Logger::error_msg("An error occurred during program execution:\n{}", e.what());
        return -1;
    }
    common::Logger::basic_msg("The calculations have been completed successfully.");
    return 0;
}