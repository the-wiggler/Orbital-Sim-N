cmake_minimum_required(VERSION 3.20)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake")

# Set vcpkg target triplet if not already set
if(NOT DEFINED VCPKG_TARGET_TRIPLET)
    set(VCPKG_TARGET_TRIPLET "x64-windows")
endif()

project(OrbitSimulation C)

# Find SDL3 packages
find_package(SDL3 CONFIG REQUIRED)
find_package(SDL3_ttf CONFIG REQUIRED)

# Source files
set(SOURCES
    src/main.c
    src/sim_calculations.c
    src/sdl_elements.c
    src/stats_window.c
)

# Header files
set(HEADERS
    src/config.h
    src/sim_calculations.h
    src/sdl_elements.h
    src/stats_window.h
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
if(WIN32)
    set_target_properties(${PROJECT_NAME} PROPERTIES WIN32_EXECUTABLE TRUE)
endif()

# Link SDL3 libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    SDL3::SDL3
    SDL3_ttf::SDL3_ttf
)

# Copy font file to build directory
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/assets/CascadiaCode.ttf"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/CascadiaCode.ttf"
)

# Copy data files to build directory
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/data"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/"
)

# On Windows, copy SDL3 DLLs to build directory
if(WIN32)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL3::SDL3>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
    )
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL3_ttf::SDL3_ttf>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
    )
endif()
