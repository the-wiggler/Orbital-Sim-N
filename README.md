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
  - Edit the JSON file `planet_data.json` to import bodies:
    ```json
    {
      "bodies": [
        {
          "name": "Earth",
          "mass": 5.972e24,
          "pos_x": 0,
          "pos_y": 0,
          "vel_x": 0,
          "vel_y": 0
        }
      ]
    }
    ```
    - **name:** the name of the body
    - **mass:** the mass (kilograms)
    - **pos_x:** the relative position in the x direction from center (meters)
    - **pos_y:** the relative position in the y direction from center (meters)
    - **vel_x:** the velocity in the x direction (meters per second)
    - **vel_y:** the velocity in the y direction (meters per second)

- **Import spacecraft**
  - Edit the JSON file `spacecraft_data.json` to import spacecraft:
    ```json
    {
      "spacecraft": [
        {
          "name": "Apollo Test Craft",
          "pos_x": -1e8,
          "pos_y": 0,
          "vel_x": 0,
          "vel_y": -1500,
          "dry_mass": 5000,
          "fuel_mass": 10000,
          "thrust": 50000,
          "specific_impulse": 300,
          "mass_flow_rate": 0.1,
          "burn_start_time": 28800,
          "burn_duration": 360,
          "burn_heading": 4.712389,
          "burn_throttle": 1.0
        }
      ]
    }
    ```

- *NOTE: these JSON files should be stored in the same directory as the program executable for the data to be read properly!

### Controls
- **Space**: Pause/resume simulation
- **R**: Reset simulation to initial state
- **Speed Control Box**: Hoover over and scroll to adjust simulation speed

### Compilation Dependencies
- **SDL3**: Graphics rendering and window management
- **SDL3_ttf**: TrueType font rendering
- **cJSON**: JSON parsing library
- **CMake**
- **vcpkg**: Package manager (used for managing dependencies)

### Font Files
- `CascadiaCode.ttf` (or modify `main.c` to use an alternative font)
