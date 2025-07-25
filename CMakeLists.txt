cmake_minimum_required(VERSION 3.17)
project(AIP)

set(CMAKE_CXX_STANDARD 14)

# set up openGL
find_package(OpenGL REQUIRED)

include_directories("include")
link_directories(lib)

# Imgui - use imgui docking branch
include_directories("Include/imgui")
include_directories("Include/imgui/backends")
set(imgui
        Include/imgui/backends/imgui_impl_glfw.cpp
        Include/imgui/backends/imgui_impl_opengl3.cpp
        Include/imgui/imgui.cpp
        Include/imgui/imgui_draw.cpp
        Include/imgui/imgui_demo.cpp
        Include/imgui/imgui_tables.cpp
        Include/imgui/imgui_widgets.cpp
        )

# STB lib
include_directories("Include/STB")
set(stb
        Include/STB/stb_image.h
        Include/STB/stb_image_write.h
        Include/STB/stb_image_resize.h
        )

# NFDE lib
add_subdirectory("Include/NFDE")

# PPMIO lib
include_directories("Include/PPMIO")
set(ppmio
        Include/PPMIO/ppm_io.h
        Include/PPMIO/ppm_io.cpp
        )


add_executable(AIP
        main.cpp
        ${imgui}
        ${stb}
        ${ppmio}
        )

target_link_libraries(AIP
        -static -static-libgcc -static-libstdc++ -static-libgcc -static-libstdc++ -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic
        OpenGL::GL
        glfw3
        glew32.lib
        glew32s.lib
        nfd
        -limm32
        )
