#pragma once

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include "glad.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "model.hpp"

struct UserPointer {
    float *zoom;
    bool *dragging;

    double *last_x;
    double *last_y;
};

class Renderer {
private:
    GLFWwindow* window;
    GLuint program;

    GLint mvp_location;
    GLint vpos_location;
    GLint vcol_location;

    GLuint vao, vbo;
public:
    float zoom = 1.0f;
    glm::vec2 cam_pos = {0.f, 0.f};

    bool dragging = false;
    double last_x = 0.0, last_y = 0.0;

    void init();
    void tick(const Frame &frame);
    void quit();

// private:
    // void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    // void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
};
