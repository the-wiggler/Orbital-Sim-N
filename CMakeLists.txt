cmake_minimum_required(VERSION 3.20)

# Project definition
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
                -mtune=native
                -flto=auto
                -funroll-loops
                -mavx2 -mfma
                -fno-math-errno -fno-trapping-math
                -DNDEBUG
        )
        add_link_options(-flto=auto)    # Enable LTO during linking

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

# macOS-specific: Add Homebrew paths for package discovery
if(APPLE)
    # Try common Homebrew installation paths
    if(EXISTS "/opt/homebrew")
        # Apple Silicon Macs
        list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")
    elseif(EXISTS "/usr/local")
        # Intel Macs
        list(APPEND CMAKE_PREFIX_PATH "/usr/local")
    endif()
endif()

# Find required packages
# Note: On Linux, ensure these are installed via package manager or built from source
find_package(SDL3 CONFIG REQUIRED)
find_package(cJSON CONFIG REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)

if(WIN32)
    find_package(PThreads4W REQUIRED)
else()
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
endif()

# Source files
set(SOURCES
        src/main.c
        src/types.h
        src/globals.h
        src/sim/bodies.h
        src/sim/bodies.c
        src/sim/spacecraft.h
        src/sim/spacecraft.c
        src/gui/SDL_engine.h
        src/gui/SDL_engine.c
        src/utility/json_loader.h
        src/utility/json_loader.c
        src/sim/simulation.h
        src/sim/simulation.c
        src/utility/telemetry_export.c
        src/utility/telemetry_export.h
        src/gui/GL_renderer.h
        src/gui/GL_renderer.c
        src/math/matrix.h
        src/gui/models.c
        src/gui/models.h
        src/gui/SDL_engine.c
        src/gui/stb_truetype.h
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Windows-specific: hide console window
if(WIN32)
    set_target_properties(${PROJECT_NAME} PROPERTIES WIN32_EXECUTABLE TRUE)
endif()

# Link libraries
# Change: Use namespaced targets (e.g., cJSON::cJSON) for Linux case-sensitivity
target_link_libraries(${PROJECT_NAME} PRIVATE
        SDL3::SDL3
        cjson
        OpenGL::GL
        GLEW::GLEW
)

if(WIN32)
    target_link_libraries(${PROJECT_NAME} PRIVATE PThreads4W::PThreads4W)
else()
    target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
    # CRITICAL FIX for Linux: Explicitly link the math library
    target_link_libraries(${PROJECT_NAME} PRIVATE m)
endif()

# Copy assets to build directory
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/data"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/"
        COMMENT "Copying data files"
)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/shaders"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/shaders"
        COMMENT "Copying shader files"
)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_SOURCE_DIR}/assets"
        "$<TARGET_FILE_DIR:${PROJECT_NAME}>/assets"
        COMMENT "Copying assets"
)

# Windows: Copy SDL3 DLLs automatically
if(WIN32)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL3::SDL3>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
            COMMENT "Copying SDL3 DLLs"
    )
endif()
