add_library(nonlinear_solver AbstractNonlinearSolver.h)

FetchContent_Declare(
        stlbfgs
        GIT_REPOSITORY https://github.com/ssloy/stlbfgs.git
        GIT_TAG "0dedf926faf66cf95fffbdc3ff5c3e0067b70296"
)
FetchContent_MakeAvailable(stlbfgs)
include_directories(${stlbfgs_SOURCE_DIR}/)
add_library(stlbfgs_adapter stlbfgs/stlbfgsAdapter.cpp stlbfgs/stlbfgsAdapter.h)
target_link_libraries(stlbfgs_adapter stlbfgs)
target_link_libraries(nonlinear_solver stlbfgs_adapter)

find_package(Eigen3 3.3)
if (Eigen3_FOUND)
    FetchContent_Declare(
            lbfgspp
            GIT_REPOSITORY https://github.com/yixuan/LBFGSpp.git
            GIT_TAG "563106b85477adcd462509d1c4c2a92cde7501f5"
    )
    FetchContent_MakeAvailable(lbfgspp)
    include_directories(${lbfgspp_SOURCE_DIR}/include/)
    add_library(lbfgspp_adapter LBFGSpp/LBFGSppAdapter.cpp LBFGSpp/LBFGSppAdapter.h)
    target_link_libraries(lbfgspp_adapter Eigen3::Eigen lbfgspp)
    target_link_libraries(nonlinear_solver lbfgspp_adapter)
endif ()

find_package(Armadillo)
if (ARMADILLO_FOUND)
    FetchContent_Declare(
            optim
            GIT_REPOSITORY https://github.com/kthohr/optim.git
            GIT_TAG "34819e0e35fcfc8ddc9028b56b9cae2c6924ae37"
    )
    FetchContent_MakeAvailable(optim)
    execute_process(COMMAND ./configure --header-only-version
            WORKING_DIRECTORY ${optim_SOURCE_DIR}/)
    include_directories(${optim_SOURCE_DIR}/header_only_version/)
    add_library(optimNM_adapter optimNM/optimNMAdapter.cpp optimNM/optimNMAdapter.h)
    target_link_libraries(optimNM_adapter ${ARMADILLO_LIBRARIES})
    target_link_libraries(nonlinear_solver optimNM_adapter)
endif ()