# rx_phys - Simple 2D Physics Engine
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

https://github.com/user-attachments/assets/a3642d7e-7ec2-4a57-aa65-fc3ef5bf20cf

**rx_phys** is a simple 2D physics demo written in a single weekend. It's designed to be lightweight, easy to understand, and serves as a demonstration of fundamental physics concepts and algorithm implementation. The core idea is to simulate the motion and interaction of rigid bodies in a 2D environment.

## Dependencies & Build Instructions

#### Dependencies:

*   C++20 compliant compiler (GCC, Clang, MSVC)
*   [sokol](https://github.com/floooh/sokol) - Graphics and utility library
*   [imgui](https://github.com/ocornut/imgui) (docking branch) - Immediate Mode GUI for visualization

#### Building:

This demo was developed using VS16 under windows, but it will build and run on linux and macos too.
To build a project it's just simple as:

```bash
$ git clone --recurse-submodules https://github.com/romixal/rx_phys
$ cd rx_phys
$ cmake -S . -B build
$ cd build
$ make
```

## Usage & Demo

#### Running the Demo:

After building, run the executable (e.g., `./rx_phys`). This will launch a window displaying a simple physics simulation. You can interact with the bodies using the mouse:

*   **Left Mouse Button:** Pick and drag bodies to move them.
*   **Right Mouse Button:** Open a context menu to add new shapes.
*   **Mouse Wheel:** Adjust zoom.
*   **Configuration Window:** Modify various simulation parameters.

## Future Work

*   Implement BVH or other spatial acceleration structure for improved performance.
*   Add support for friction and rotation.
*   Introduce kinematic bodies.
*   Explore Structure of Arrays (SOA) data layout for optimization.
*   Add spring joints.

## License

This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.
