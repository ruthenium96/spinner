#include <iostream>

#include "src/common/runner/Executer.h"
#include "src/input/Parser.h"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        std::cout << "ERROR: missing input file" << std::endl;
        return 0;
    }
    if (argc > 2) {
        std::cout << "ERROR: Too much arguments" << std::endl;
        return 0;
    }
    auto parser = input::Parser(argv[1]);
    runner::Executer::execute(parser);
    return 0;
}