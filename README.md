# Orbital Simulation Software N

A real-time gravitational physics simulator built in C with SDL3.

## Overview

This program simulates gravitational interactions between celestial bodies using Newton's law of universal gravitation.

## Features

- Real-time Physics: Gravitational force calculations using F = GMm/r^2 with configurable time steps
- CSV Data Import: Load predefined planetary systems from CSV files
- Visual UI Elements: 
  - Real-time statistics display showing body properties
  - Scale reference bar for distance measurement
  - Speed control interface
  - Pause/resume functionality

## Usage
### Simulation Setup
- **Import gravitational bodies**
  - Edit the CSV file according to the attached template to import
    - It must have filename "planet_data.csv" and must follow the template for the sim to work.
    - **Planet name:** the name of the planet you're adding - self explanatory
    - **mass:** the mass of the planet (kilograms)
    - **pos_x:** the relative position of the body in the x direction from the "center" (meters)
    - **pos_y:** the relative position of the body in the y direction from the "center" (meters)
    - **vel_x:** the "relative speed" of the body in the x direction to empty space (meters per second)
    - **vel_y:** the "relative speed" of the body in the y direction to empty space (meters per second)
  - Once this file is complete, you can run the executable file which imports data from the CSV into the simulation

### Controls
- **Space**: Pause/resume simulation
- **Left Click**: Add a new orbital body at cursor position
- **R**: Reset simulation to initial state
- **Speed Control Box**: Click and drag to adjust simulation speed

### Compilation Dependencies
- **SDL3**: Graphics rendering and window management
- **SDL3_ttf**: TrueType font rendering
- **CMake**

### Font Files
- `CascadiaCode.ttf` (or modify `main.c` to use an alternative font)
