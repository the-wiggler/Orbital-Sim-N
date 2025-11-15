cmake_minimum_required(VERSION 3.20)
project(orbit_simulation C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Build output directories - everything goes to build root
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Find SDL3 package
find_package(SDL3 REQUIRED CONFIG REQUIRED COMPONENTS SDL3)
find_package(SDL3_ttf REQUIRED CONFIG REQUIRED)

# Source files
set(SOURCES
    src/main.c
    src/calculation_functions.c
    src/sdl_elements.c
)

# Header files (for IDE organization)
set(HEADERS
    src/config.h
    src/calculation_functions.h
    src/sdl_elements.h
)

# Create executable
add_executable(orbital_sim ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(orbital_sim
    PRIVATE
    SDL3::SDL3
    SDL3_ttf::SDL3_ttf
)

# Include directories
target_include_directories(orbital_sim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Copy font file directly to build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/assets/CascadiaCode.ttf
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CascadiaCode.ttf COPYONLY)

# Copy CSV file directly to build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/planet_data.csv
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/planet_data.csv COPYONLY)

# Platform-specific settings
if(WIN32)
    # Windows-specific settings
    target_link_libraries(orbital_sim PRIVATE m)
elseif(UNIX)
    # Unix/Linux-specific settings
    target_link_libraries(orbital_sim PRIVATE m)
endif()

# Enable warnings
if(MSVC)
    target_compile_options(orbital_sim PRIVATE /W4)
else()
    target_compile_options(orbital_sim PRIVATE -Wall -Wextra -Wpedantic)
endif()
