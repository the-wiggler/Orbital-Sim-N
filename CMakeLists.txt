cmake_minimum_required(VERSION 3.20)

# Project definition
project(OrbitSimulation C)

# C Standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- 1. FIND PACKAGES ---
find_package(SDL3 CONFIG REQUIRED)
find_package(cJSON CONFIG REQUIRED)

if(NOT EMSCRIPTEN)
    find_package(OpenGL REQUIRED)
    find_package(GLEW REQUIRED)
endif()

if(WIN32)
    find_package(PThreads4W REQUIRED)
else()
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
if(WIN32)
    set_target_properties(${PROJECT_NAME} PROPERTIES WIN32_EXECUTABLE TRUE)
endif()

# --- 4. COMPILER FLAGS & OPTIMIZATIONS ---
if(EMSCRIPTEN)
    # --- WEB ASSEMBLY FLAGS ---
    target_compile_options(${PROJECT_NAME} PRIVATE -msimd128)

    # Linker flags for WebGL 2
    target_link_options(${PROJECT_NAME} PRIVATE
            "-sUSE_WEBGL2=1"
            "-sFULL_ES3=1"
            "-sALLOW_MEMORY_GROWTH=1" # Recommended for sims
    )

    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".html")

elseif(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    # --- DESKTOP (GCC/CLANG) FLAGS ---
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)

    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(${PROJECT_NAME} PRIVATE
                -O3
                -march=native
                -mtune=native
                -flto=auto
                -funroll-loops
                -mavx2 -mfma
                -fno-math-errno -fno-trapping-math
                -DNDEBUG
        )
        target_link_options(${PROJECT_NAME} PRIVATE -flto=auto)

        if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
            target_compile_options(${PROJECT_NAME} PRIVATE -fno-stack-protector)
        endif()
    endif()

    if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        target_compile_options(${PROJECT_NAME} PRIVATE -O1 -g -fno-inline)
    endif()

elseif(MSVC)
    # --- WINDOWS (MSVC) FLAGS ---
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(${PROJECT_NAME} PRIVATE
                /O2
                /GL
                /arch:AVX2
                /fp:fast
        )
        target_link_options(${PROJECT_NAME} PRIVATE /LTCG)
    endif()
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

if(WIN32)
    target_link_libraries(${PROJECT_NAME} PRIVATE PThreads4W::PThreads4W)
else()
    target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads m)
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

if(WIN32)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL3::SDL3>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
            COMMENT "Copying SDL3 DLLs"
    )
endif()