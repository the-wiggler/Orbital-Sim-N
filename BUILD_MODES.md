# Orbit Simulation N - Build Modes

This project supports two build modes: **GUI mode** (default) and **Headless mode** (command-line).

## GUI Mode (Default)

Build with graphical interface, OpenGL rendering, and SDL3:

```bash
cmake -B build
cmake --build build
./build/OrbitSimulation
```

## Headless Mode (Command-Line)

Build without GUI dependencies for terminal/server use:

```bash
cmake -B build -DENABLE_GUI=OFF
cmake --build build
./build/OrbitSimulation
```

### Headless Commands

When running in headless mode, you can control the simulation with these commands:

- `start` or `run` - Start the simulation
- `stop` or `pause` - Pause the simulation
- `status` - Show current simulation state
- `reset` - Reset the simulation
- `load` - Load simulation parameters from JSON
- `step <num>` - Change simulation time step (e.g., `step 0.01`)
- `logfreq <num>` - Change CSV logging frequency (e.g., `logfreq 0.1`)
- `quit` or `exit` - Exit the program
- `help` - Show available commands

## File Structure

- `src/main.c` - GUI main function and shared physics simulation thread
- `src/headless_main.c` - Command-line interface and headless main function
- Build system automatically selects the appropriate files based on `ENABLE_GUI` option
