add_subdirectory(common)
add_subdirectory(entities)
add_subdirectory(group)
add_subdirectory(model)
add_subdirectory(nonlinear_solver)
add_subdirectory(space)
add_subdirectory(spin_algebra)

add_executable(main_components
        main/main_components.cpp
        )

target_link_libraries(main_components space group space common)
