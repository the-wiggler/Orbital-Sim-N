cmake_minimum_required(VERSION 3.20)

# Project definition (should come before compiler settings)
project(OrbitSimulation C)

# C Standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler-specific flags
if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        # Aggressive optimization flags for maximum performance
        add_compile_options(
            -O3
            -march=native
            -flto
            -ffast-math
            -funroll-loops
        )
        add_link_options(-flto)    # Enable LTO during linking

        # GCC-specific optimizations
        if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
            add_compile_options(-fno-stack-protector)  # Remove stack protection overhead
        endif()
    endif()
    # For profiling: keep debug symbols and prevent aggressive inlining
    if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        add_compile_options(-O1 -g -fno-inline)
    endif()
elseif(MSVC)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        # MSVC optimization flags
        add_compile_options(
            /O2                    # Maximum optimization
            /GL                    # Whole program optimization
            /arch:AVX2             # Use AVX2 instructions if available
            /fp:fast               # Fast floating-point model
        )
        add_link_options(/LTCG)    # Link-time code generation
    endif()
endif()

# Find required packages
find_package(SDL3 CONFIG REQUIRED)
find_package(SDL3_ttf CONFIG REQUIRED)
find_package(cJSON CONFIG REQUIRED)

if(WIN32)
    find_package(PThreads4W REQUIRED)
else()
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
endif()

# Source files (no need for separate HEADERS variable)
set(SOURCES
        src/main.c
        src/types.h
        src/globals.h
        src/globals.c
        src/bodies.h
        src/bodies.c
        src/spacecraft.h
        src/spacecraft.c
        src/renderer.h
        src/renderer.c
        src/json_loader.h
        src/json_loader.c
        src/simulation.h
        src/simulation.c
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Windows-specific: hide console window
if(WIN32)
    set_target_properties(${PROJECT_NAME} PROPERTIES WIN32_EXECUTABLE TRUE)
endif()

# Link libraries (unified approach)
target_link_libraries(${PROJECT_NAME} PRIVATE
        SDL3::SDL3
        SDL3_ttf::SDL3_ttf
        cjson
)

if(WIN32)
    target_link_libraries(${PROJECT_NAME} PRIVATE PThreads4W::PThreads4W)
else()
    target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
endif()

# Copy assets to build directory
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/assets/font.ttf"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/font.ttf"
        COMMENT "Copying font file"
)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/data"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/"
        COMMENT "Copying data files"
)

# Windows: Copy SDL3 DLLs automatically
if(WIN32)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL3::SDL3>
            $<TARGET_FILE:SDL3_ttf::SDL3_ttf>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
            COMMENT "Copying SDL3 DLLs"
    )
endif()