cmake_minimum_required(VERSION 3.16)
project(app CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(app
    app.cpp
    glad/src/glad.c
    src/main.cpp
    src/utils.cpp
    src/model/model.cpp
    src/model/neighborhood_search/neighborhood_search.cpp
    src/model/physics/elastic_solids.cpp
    src/model/physics/gravity.cpp
    src/model/physics/liquid.cpp
    src/renderer/renderer.cpp
)

option(ENABLE_ASAN_UBSAN "Build with AddressSanitizer + UndefinedBehaviorSanitizer instead of the optimized build" OFF)

if(ENABLE_ASAN_UBSAN)
    set(CMAKE_BUILD_TYPE Debug)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)

    target_compile_options(app PRIVATE
        -g -Og -fno-omit-frame-pointer
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wcast-align -Wunused -Wnull-dereference -Wdouble-promotion
        -Wformat=2 -Wimplicit-fallthrough -Wnon-virtual-dtor -Wold-style-cast
        -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast
        -fsanitize=address,undefined,leak
        -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow
        -fsanitize=integer-divide-by-zero -fsanitize=null -fsanitize=alignment
        -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=object-size
        -fsanitize=vptr -fsanitize=nonnull-attribute -fsanitize=returns-nonnull-attribute
        -fsanitize=unreachable -fsanitize=signed-integer-overflow -fsanitize=shift
        -fsanitize=pointer-overflow -fno-sanitize-recover=all
    )
    target_link_options(app PRIVATE -fsanitize=address,undefined,leak)
else()
    set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

    target_compile_options(app PRIVATE
        -O3
        -march=native
        -mtune=native
        -fopenmp
        -ffast-math
        -fno-math-errno
        -funsafe-math-optimizations
        -funroll-loops
        -fomit-frame-pointer
        -Wall
        -Wextra
    )
    target_link_options(app PRIVATE -flto)
endif()

target_link_libraries(app PRIVATE OpenMP::OpenMP_CXX)

target_include_directories(app PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/glad/include
    ${CMAKE_SOURCE_DIR}/glad/include/glad
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/model
    ${CMAKE_SOURCE_DIR}/src/model/physics
    ${CMAKE_SOURCE_DIR}/src/model/neighborhood_search
    ${CMAKE_SOURCE_DIR}/src/model/barnes_hut
    ${CMAKE_SOURCE_DIR}/src/renderer
)

target_link_libraries(app PRIVATE
    OpenGL::GL
    glfw
    glm::glm
)