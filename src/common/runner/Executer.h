#ifndef SPINNER_EXECUTER_H
#define SPINNER_EXECUTER_H

#include "Runner.h"
#include "src/input/Parser.h"

namespace spinner::runner {

class Executer {
  public:
    static void execute(input::Parser parser);
    static void dry_execute(input::Parser parser);
};

} // namespace spinner::runner

#endif  //SPINNER_EXECUTER_H
