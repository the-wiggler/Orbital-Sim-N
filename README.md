[![Nightly Native Builds](https://github.com/toastyy-1/Orbital-Sim-N/actions/workflows/matrix-deploy.yml/badge.svg?branch=main)](https://github.com/toastyy-1/Orbital-Sim-N/actions/workflows/matrix-deploy.yml) [![WASM Build](https://github.com/toastyy-1/Orbital-Sim-N/actions/workflows/web-deploy.yml/badge.svg?branch=main)](https://github.com/toastyy-1/Orbital-Sim-N/actions/workflows/web-deploy.yml)

# Orbital Simulation Software N

A real-time 3D gravitational physics simulator built in C with OpenGL and SDL3.

### Gravitational Physics
- Real-time gravitational force calculations using F = GMm/r²
- Verlet Integration
- Collision Detection
- Adjustable simulation speed

### Spacecraft Systems
- **Propulsion**:
  - Thrust with configurable specific impulse
  - Fuel consumption based on mass flow rate
  - Variable throttle control (0-100%)

- **Attitude Control**:
  - Quaternion-based 3D orientation system
  - Automatic attitude control during burns
  - Configurable moment of inertia and nozzle gimbal range

- **Burn Planning**:
  - Schedule multiple burns with precise start times and durations
    - **Tangent**: Aligns thrust with velocity vector (prograde/retrograde burns)
    - **Absolute**: Burns relative to the inertial reference frame (space coordinates)
    - **Normal**: Perpendicular to orbit -- *currently unsupported :(*
  - Each burn type can be configured with a rotation heading offset

## Usage

### Controls

#### Mouse Controls
- **Right Mouse Button + Drag**: Rotate camera around the scene
- **Mouse Wheel Up**: Zoom in
- **Mouse Wheel Down**: Zoom out

#### Console Commands
The program features a command-line interface at the bottom of the window. Available commands:

| Command | Description |
|---------|-------------|
| `load` | Load simulation from `simulation_data.json` |
| `pause` or `p` | Pause the simulation |
| `resume` or `r` | Resume the simulation |
| `reset` | Reset the simulation to initial state |
| `step <value>` | Set simulation time step (e.g., `step 0.01`) |
| `enable guidance-lines` | Show lines between celestial bodies |
| `disable guidance-lines` | Hide lines between celestial bodies |

**Note**: Type commands in the console at the bottom of the window and press Enter to execute.

### Configuration Files

The simulation is configured via `simulation_data.json`:

#### Adding Celestial Bodies

```json
{
  "bodies": [
    {
      "name": "Earth",
      "mass": 5.972e24,
      "pos_x": 0,
      "pos_y": 0,
      "pos_z": 0,
      "vel_x": 0,
      "vel_y": 0,
      "vel_z": 0,
      "radius": 6371000
    }
  ]
}
```

**Parameters:**
- `name`: Display name
- `mass`: Mass in kilograms
- `pos_x`, `pos_y`, `pos_z`: Position in meters (m) from the origin
- `vel_x`, `vel_y`, `vel_z`: Initial velocity in meters per second
- `radius`: Body radius in meters

#### Adding Spacecraft

```json
{
  "spacecraft": [
    {
      "name": "Example Craft",
      "pos_x": 6671000.0,
      "pos_y": 0.0,
      "pos_z": 0.0,
      "vel_x": 0.0,
      "vel_y": 7700.0,
      "vel_z": 0.0,
      "dry_mass": 5500.0,
      "fuel_mass": 18000.0,
      "thrust": 450000.0,
      "specific_impulse": 314.0,
      "mass_flow_rate": 146.0,
      "attitude": 0.0,
      "moment_of_inertia": 20000.0,
      "nozzle_gimbal_range": 0.105,
      "burns": [
        {
          "burn_target": "Earth",
          "burn_type": "tangent",
          "start_time": 3600.0,
          "duration": 1000.0,
          "heading": 0.0,
          "throttle": 1.0
        }
      ]
    }
  ]
}
```

**Spacecraft Parameters:**
- **Position/Velocity**: 3D coordinates (m) and velocities (m/s)
  - `pos_x`, `pos_y`, `pos_z`: Initial position
  - `vel_x`, `vel_y`, `vel_z`: Initial velocity
- **Mass Properties**:
  - `dry_mass`: Spacecraft mass without fuel (kg)
  - `fuel_mass`: Available propellant (kg)
- **Propulsion**:
  - `thrust`: Maximum engine thrust (N)
  - `specific_impulse`: Engine efficiency (s)
  - `mass_flow_rate`: Fuel consumption rate (kg/s)
- **Attitude**:
  - `attitude`: Initial rotation angle about Z-axis (radians)
  - `moment_of_inertia`: Rotational inertia (kg⋅m²)
  - `nozzle_gimbal_range`: Thrust vectoring limit (radians)

**Burn Parameters:**
- `burn_target`: Name of reference body (must match body name exactly -- case sensitive)
- `burn_type`: Burn reference frame type
  - `"tangent"`: Thrust aligned with velocity vector relative to target body
  - `"absolute"`: Thrust in absolute inertial reference frame
  - `"normal"`: *Currently unsupported*
- `start_time`: Burn start time (seconds from T+0)
- `duration`: Burn length (seconds)
- `heading`: Rotation offset from base burn direction (radians)
- `throttle`: Engine power (0.0 to 1.0)

Note: Due to the nature of time step based simulation, the craft burn time may not start on the "exact" time that it is expected to. Variations in accuracy in this regard can be attributed having too high of a time step during burns.

## Building

### Dependencies
- **SDL3**: Window management and OpenGL context creation
- **OpenGL**: 3D graphics rendering
- **GLEW**: OpenGL extension loading library
- **cJSON**: JSON parsing library
- **pthreads**: Multi-threading support (or PThreads4W on Windows)
  - Note: better Windows multi-threading support will be added in the future. Major slowdowns on Windows should be expected.

### Required Files
- **Shader files**: `shaders/simple.vert` and `shaders/simple.frag`
- **Font file**: `font.ttf` (placed in build directory)
- **Configuration**: `simulation_data.json`

### Build Instructions
Use CMake to build the project:

The build system automatically copies required assets (shaders, fonts, data files) to the build directory.

### Web Build Instructions
#### Build with Conan Dependencies for Web
```sh
mkdir build
cd build
conan install .. --build=missing -s build_type=Release -pr ../web/conan_profile
# Use either
# Windows:
#  call .\Release\generators\conanbuild.bat
#  cmake .. -G "MinGW Makefiles" -DCMAKE_TOOLCHAIN_FILE=Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
# Or
# Linux:
#  source ./Release/generators/conanbuild.sh
#  cmake .. -DCMAKE_TOOLCHAIN_FILE=Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

#### Quick local server with correct headers
For local testing of the WebAssembly build with pthreads, the repository includes a small helper server that sets the required COOP/COEP headers and opens the page automatically.

Usage examples:
- `python ./web/serve_web.py` — auto-detects the build folder and opens OrbitSimulation.html on 127.0.0.1:8000
- `python ./web/serve_web.py --port 8080` — choose a different port

The server adds these headers on every response so that threads are enabled in the browser:
- Cross-Origin-Opener-Policy: same-origin
- Cross-Origin-Embedder-Policy: require-corp
