#pragma once

#include "entities/Space.h"
#include <memory>

using TaskPtr = std::shared_ptr<Space>;

class Processor{
    virtual TaskPtr process(TaskPtr task) = 0;
};