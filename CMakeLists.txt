cmake_minimum_required(VERSION 3.17)
project(Cusa)

##===---------------------------------------------------------------------===##
## Sokol Wrapper
##===---------------------------------------------------------------------===##

# Sokol
add_library(sokol STATIC "lib/sokol_wrapper/sokol.cpp")
if(EMSCRIPTEN)
    target_link_libraries(sokol PUBLIC sokol_headers)
elseif(APPLE) # Apple specific libraries
    target_link_libraries(sokol
        PUBLIC
            sokol_headers
        PRIVATE
            "-framework AudioToolbox"
            "-framework Cocoa"
            "-framework QuartzCore"
            "-framework Metal"
            "-framework MetalKit"
    )
    target_compile_options(sokol PRIVATE -x objective-c)
elseif(UNIX) # Unix-specific libraries
    find_package(X11 REQUIRED)
    find_package(OpenGL REQUIRED)
    target_link_libraries(sokol
        PUBLIC
            sokol_headers
        PRIVATE
            OpenGL::GL
            X11
            Xi
            Xcursor
    )
else()
    target_link_libraries(sokol PUBLIC sokol_headers)
endif()
target_include_directories(sokol PUBLIC "lib/sokol_wrapper/")

##===---------------------------------------------------------------------===##
## Main App
##===---------------------------------------------------------------------===##

add_executable(rx_phys rx_phys.cpp)
# NOTE(kobylinsky): C++20 mainly used for designated initializers.
set_property(TARGET rx_phys PROPERTY CXX_STANDARD 20)
target_link_libraries(rx_phys PRIVATE sokol imgui)

##===---------------------------------------------------------------------===##
## External Libraries
##===---------------------------------------------------------------------===##

# Sokol libraries
add_library(sokol_headers INTERFACE)
target_include_directories(sokol_headers INTERFACE "external/sokol")

# Dear ImGui
add_library(imgui STATIC
    "external/imgui/imgui.cpp"
    "external/imgui/imgui_draw.cpp"
    "external/imgui/imgui_widgets.cpp"
    "external/imgui/imgui_tables.cpp"
    "external/imgui/imgui_demo.cpp"
    "external/imgui/imgui.h"
    "external/imgui/imgui_internal.h"
    "external/imgui/misc/cpp/imgui_stdlib.cpp"
)
target_include_directories(imgui PUBLIC "external/imgui")
target_compile_definitions(imgui PUBLIC
    IMGUI_DEFINE_MATH_OPERATORS
    IMGUI_DISABLE_OBSOLETE_FUNCTIONS
)