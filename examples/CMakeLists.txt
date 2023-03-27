add_executable(10.1021_ja0305959
        10.1021_ja0305959/d-nit2.h
        10.1021_ja0305959/main.cpp)

target_link_libraries(10.1021_ja0305959 common)

add_executable(10.1039_B414284E
        10.1039_B414284E/3d-nit2.h
        10.1039_B414284E/dtd-nit2.h
        10.1039_B414284E/main.cpp)

target_link_libraries(10.1039_B414284E common)

add_executable(10.1039_D0DT03184D
        10.1039_D0DT03184D/manganese_cocrystal.h
        10.1039_D0DT03184D/nickel_crystal.h
        10.1039_D0DT03184D/main.cpp)

target_link_libraries(10.1039_D0DT03184D common)
