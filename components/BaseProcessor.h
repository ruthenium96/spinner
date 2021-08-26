#pragma once

#include "common/Task.h"
#include <memory>

using TaskPtr = std::shared_ptr<Task>;

class BaseProcessor{
    virtual TaskPtr process(TaskPtr task) = 0;
};