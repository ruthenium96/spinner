add_library(common
        LexicographicIndexConverter.cpp LexicographicIndexConverter.h
        Logger.h
        Runner.cpp)

target_link_libraries(common components)