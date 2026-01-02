# Orbital Simulation Software N

A real-time 3D gravitational physics simulator built in C with OpenGL and SDL3.

### Gravitational Physics
- Real-time gravitational force calculations using F = GMm/r²
- Verlet Integration
- Total system energy tracking with error monitoring
- Collision Detection
- Adjustable simulation speed

!!!!NEEDS TO BE UPDATED:
###  Spacecraft Systems
  - Thrust  with configurable specific impulse
  - Fuel consumption based on mass flow rate
  - Variable throttle control (0-100%)

- **Burn Planning**:
  - Schedule multiple burns with precise start times and durations
    - **Tangent**: retrograde or prograde burns
    - **Normal**: burns perpendicular to the orbit
    - **Absolute**: burns relative to the inertial observer (space)
  - Each of these burn types can be configured with a rotation heading relative to each respective axis

## Usage
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
- `pos_x`, `pos_y`: Position in meters (m) from the "center" (center of the screen)
- `vel_x`, `vel_y`: Initial velocity in meters per second
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
- **Position/Velocity**: Same as bodies (m, m/s)
- **Mass Properties**:
  - `dry_mass`: Spacecraft mass without fuel (kg)
  - `fuel_mass`: Available propellant (kg)
- **Propulsion**:
  - `thrust`: Maximum engine thrust (N)
  - `specific_impulse`: Engine efficiency (s)
  - `mass_flow_rate`: Fuel consumption rate (kg/s)
- **Attitude**:
  - `attitude`: Initial orientation (radians)
  - `moment_of_inertia`: Rotational inertia (kg⋅m^2)
  - `nozzle_gimbal_range`: Thrust vectoring limit (radians)

**Burn Parameters:**
- `burn_target`: Name of reference body (must match body name exactly -- case sensitive)
- `burn_type`: `"tangent"`, `"normal"`, or `"absolute"`
- `start_time`: Burn start time (seconds from T+0)
- `duration`: Burn length (seconds)
- `heading`: Direction offset (radians)
- `throttle`: Engine power (0.0 to 1.0)

## Building

### Dependencies
- **SDL3**: Window management and OpenGL context creation
- **OpenGL**: 3D graphics rendering
- **GLEW**: OpenGL extension loading library
- **FreeType**: TrueType font rendering for on-screen text
- **cJSON**: JSON parsing library
- **pthreads**: Multi-threading support (or PThreads4W on Windows)

### Required Files
- **Shader files**: `shaders/simple.vert` and `shaders/simple.frag`
- **Font file**: `font.ttf` (placed in build directory)
- **Configuration**: `simulation_data.json`

### Build Instructions
Use CMake to build the project:

The build system automatically copies required assets (shaders, fonts, data files) to the build directory.

### Build with Conan Dependencies for Web
```sh
mkdir build && cd build
conan install .. --build=missing -s build_type=Release -pr ../web
call .\Release\generators\conanbuild.bat
cmake .. -G "MinGW Makefiles" -DCMAKE_TOOLCHAIN_FILE=Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

### Quick local server with correct headers
For local testing of the WebAssembly build with pthreads, the repository includes a small helper server that sets the required COOP/COEP headers and opens the page automatically.

Usage examples:
- `python serve_web.py` — auto-detects the build folder and opens OrbitSimulation.html on 127.0.0.1:8000
- `python serve_web.py --root build/Release` — specify the directory that contains OrbitSimulation.html
- `python serve_web.py --port 8080` — choose a different port

The server adds these headers on every response so that threads are enabled in the browser:
- Cross-Origin-Opener-Policy: same-origin
- Cross-Origin-Embedder-Policy: require-corp
