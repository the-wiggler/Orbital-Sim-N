cmake_minimum_required(VERSION 3.20)

# Project definition
project(OrbitSimulation C)

# C Standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Force Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Find required packages
# Note: On Linux, ensure these are installed via package manager or built from source
# --- 1. FIND PACKAGES ---
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

find_package(SDL3 CONFIG REQUIRED)
find_package(cJSON CONFIG REQUIRED)

if(NOT EMSCRIPTEN)
    find_package(OpenGL REQUIRED)
    find_package(GLEW REQUIRED)
endif()

if(NOT EMSCRIPTEN)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
endif()

# --- 2. DEFINE SOURCES ---
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
        src/gui/stb_truetype.h
)

# --- 3. CREATE EXECUTABLE ---
add_executable(${PROJECT_NAME} ${SOURCES})

# Windows-specific: hide console window
#if(WIN32)
#    set_target_properties(${PROJECT_NAME} PROPERTIES WIN32_EXECUTABLE TRUE)
#endif()

# --- 4. COMPILER FLAGS & OPTIMIZATIONS ---
if(EMSCRIPTEN)
    # --- WEB ASSEMBLY FLAGS ---
    target_compile_options(${PROJECT_NAME} PRIVATE -msimd128 -pthread)

    # Linker flags for WebGL 2
    target_link_options(${PROJECT_NAME} PRIVATE
            "-sUSE_WEBGL2=1"
            "-sFULL_ES3=1"
            "-sUSE_PTHREADS=1"
            "-sASYNCIFY=1"
            "-pthread"

            # Mappings for virtual filesystem
            "SHELL:--preload-file ${CMAKE_SOURCE_DIR}/assets@/assets"
            "SHELL:--preload-file ${CMAKE_SOURCE_DIR}/shaders/web@/shaders"
            "SHELL:--preload-file ${CMAKE_SOURCE_DIR}/data@/"
    )

    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".html")

elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-function)
    # Aggressive optimization flags for maximum performance
    add_compile_options(
            -O3
            -march=native
            -mtune=native
            -flto=auto
            -funroll-loops
            -fno-math-errno -fno-trapping-math
            -DNDEBUG
    )

    # Architecture-specific SIMD optimizations
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i686")
        # x86/x64: Use AVX2 and FMA instructions
        add_compile_options(-mavx2 -mfma)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64|ARM64")
        # ARM: Use NEON instructions (enabled by -march=native on ARM)
        # No additional flags needed as -march=native handles this
    endif()

    add_link_options(-flto=auto)    # Enable LTO during linking

    # GCC-specific optimizations
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
        add_compile_options(-fno-stack-protector)  # Remove stack protection overhead
    endif()
elseif(MSVC)
    # MSVC optimization flags (matched to Linux -O3 level)
    add_compile_options(
            /W4                    # High warning level
            /permissive-           # Standards conformance mode (pedantic)
            /Ox                    # Maximum optimization (aggressive, similar to -O3)
            /Ot                    # Favor fast code
            /Oi                    # Generate intrinsic functions
            /GL                    # Whole program optimization
            /arch:AVX2             # Use AVX2 instructions if available
            /fp:fast               # Fast floating-point model
            /GS-                   # Disable security checks for performance
    )
    add_link_options(/LTCG)    # Link-time code generation
endif()

# --- 5. INCLUDE PATHS (macOS) ---
if(APPLE)
    if(EXISTS "/opt/homebrew")
        target_include_directories(${PROJECT_NAME} PRIVATE "/opt/homebrew/include")
        link_directories("/opt/homebrew/lib")
    elseif(EXISTS "/usr/local")
        target_include_directories(${PROJECT_NAME} PRIVATE "/usr/local/include")
        link_directories("/usr/local/lib")
    endif()
endif()

# --- 6. LINK LIBRARIES ---
set(EXTRA_LIBS "")
if(NOT EMSCRIPTEN)
    list(APPEND EXTRA_LIBS OpenGL::GL GLEW::GLEW)
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE
        SDL3::SDL3
        cjson
        ${EXTRA_LIBS}
)

if(NOT EMSCRIPTEN)
    if(WIN32)
        target_link_libraries(${PROJECT_NAME} PRIVATE PThreads4W::PThreads4W)
    else()
        target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
        # CRITICAL FIX for Linux: Explicitly link the math library
        target_link_libraries(${PROJECT_NAME} PRIVATE m)
    endif()
endif()

# --- 7. ASSET COPYING ---
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
