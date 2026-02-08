#include <gtest/gtest.h>

#include "src/common/Logger.h"

int main(int argc, char** argv) {
    spinner::common::Logger::set_level(spinner::common::PrintLevel::off);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}