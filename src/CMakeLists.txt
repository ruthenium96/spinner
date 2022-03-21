add_subdirectory(group)
add_subdirectory(common)
add_subdirectory(entities)
add_subdirectory(components)

add_executable(main_components
        main/main_components.cpp
        )

target_link_libraries(main_components components group space)
