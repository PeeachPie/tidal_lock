cmake_minimum_required(VERSION 3.16)
project(app CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)

add_executable(app
    app.cpp
    glad/src/glad.c
    src/main.cpp
    src/utils.cpp
    src/model/model.cpp
    src/renderer/renderer.cpp
)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

target_compile_options(app PRIVATE
    -O3
    -march=native
    -mtune=native
    -ffast-math
    -fno-math-errno
    -funsafe-math-optimizations
    -funroll-loops
    -fomit-frame-pointer
)

target_link_options(app PRIVATE
    -flto
)

find_package(OpenMP REQUIRED)

target_link_libraries(app PRIVATE OpenMP::OpenMP_CXX)

target_include_directories(app PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/glad/include
    ${CMAKE_SOURCE_DIR}/glad/include/glad
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/model
    ${CMAKE_SOURCE_DIR}/src/renderer
    ${CMAKE_SOURCE_DIR}/src/barnes_hut
)

target_link_libraries(app PRIVATE
    OpenGL::GL
    glfw
    glm::glm
)