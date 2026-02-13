#include "src/common/Logger.h"
#include "src/common/runner/Executer.h"
#include "src/input/Parser.h"

const std::string help_string = "\nUsage: ./spinner_main filename.yml [--dry]";

int main(int argc, char *argv[]) {
    spinner::common::Logger::set_pattern("[%C/%m/%d %H:%M:%S.%e] [%t] %v");
    if (argc == 1) {
        spinner::common::Logger::error_msg("Missing input file.{}", help_string);
        return -1;
    }
    if (argc > 3) {
        spinner::common::Logger::error_msg("Too much arguments.{}", help_string);
        return -1;
    }
    bool dry_run = false;
    if (argc == 3 && argv[2] != std::string("--dry")) {
        spinner::common::Logger::error_msg("Unknown argument: {}{}", argv[2], help_string);
        return -1;
    }
    if (argc == 3 && argv[2] == std::string("--dry")) {
        dry_run = true;
    }
    try {
        if (dry_run) {
            auto parser = spinner::input::Parser(argv[1], dry_run);
            spinner::runner::Executer::dry_execute(parser);
        } else {
            spinner::common::Logger::basic_msg("The calculations have started");
            spinner::common::Logger::separate(0, spinner::common::basic);        
            auto parser = spinner::input::Parser(argv[1]);
            spinner::runner::Executer::execute(parser);
        }
    } catch (const std::exception& e) {
        spinner::common::Logger::error_msg("An error occurred during program execution:\n{}", e.what());
        return -1;
    }
    spinner::common::Logger::basic_msg("The calculations have been completed successfully.");
    return 0;
}