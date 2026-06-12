#include "renderer.hpp"

struct Vertex {
    glm::vec2 pos;
    glm::vec3 col;
};

const char* vertex_shader_text =
    "#version 330 core\n"
    "uniform mat4 MVP;\n"
    "in vec2 vPos;\n"
    "in vec3 vCol;\n"
    "out vec3 color;\n"
    "void main() {\n"
    "    gl_Position = MVP * vec4(vPos, 0.0, 1.0);\n"
    "    color = vCol;\n"
    "}\n";

const char* fragment_shader_text =
    "#version 330 core\n"
    "in vec3 color;\n"
    "out vec4 fragColor;\n"
    "void main() {\n"
    "    fragColor = vec4(color, 1.0);\n"
    "}\n";

void appendCircle(std::vector<Vertex>& out,
                  float cx, float cy, float r,
                  glm::vec3 color,
                  int& start, int& count,
                  int segments = 32)
{
    start = out.size();

    out.push_back({{cx, cy}, color});

    for (int i = 0; i <= segments; i++) {
        float a = 2.f * M_PI * i / segments;

        float x = cx + r * cos(a);
        float y = cy + r * sin(a);

        out.push_back({{x, y}, color});
    }

    count = (segments + 2);
}

void appendPoints(std::vector<Vertex>& out, const std::vector<ParticleFrame>& particles, int& start, int& count)
{
    start = out.size();

    for (auto& p : particles) {
        out.push_back({
            {(float)p.pos.x, (float)p.pos.y},
            {1.f, 1.f, 1.f}
        });
    }

    count = particles.size();
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS || action == GLFW_REPEAT)
    {
        if (mods & GLFW_MOD_CONTROL)
        {
            float& zoom = ((Renderer*)glfwGetWindowUserPointer(window))->zoom;

            if (key == GLFW_KEY_EQUAL)   // Ctrl +
                zoom *= 0.9f;

            if (key == GLFW_KEY_MINUS)   // Ctrl -
                zoom *= 1.f / 0.9f;
        }
    }

    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        Renderer* r = (Renderer*)glfwGetWindowUserPointer(window);
        
        if (action == GLFW_PRESS)
        {
            r->dragging = true;
            glfwGetCursorPos(window, &r->last_x, &r->last_y);
        }
        else if (action == GLFW_RELEASE)
        {
            r->dragging = false;
        }
    }
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    Renderer* r = (Renderer*)glfwGetWindowUserPointer(window);

    if (!r->dragging) return;

    double dx = xpos - r->last_x;
    double dy = ypos - r->last_y;

    r->last_x = xpos;
    r->last_y = ypos;

    int width, height;
    glfwGetWindowSize(window, &width, &height);

    float ratio = width / (float)height;

    float world_dx = -dx / width  * 2.0f * ratio * r->zoom;
    float world_dy =  dy / height * 2.0f * r->zoom;

    r->cam_pos.x += world_dx;
    r->cam_pos.y += world_dy;
}

void Renderer::init() {
    glfwInit();

    window = glfwCreateWindow(1600, 1200, "Приливы", NULL, NULL);
    glfwMakeContextCurrent(window);
    gladLoadGL();

    glfwSetWindowUserPointer(window, this);
    glfwSetKeyCallback(window, key_callback);

    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertex_shader_text, NULL);
    glCompileShader(vs);

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragment_shader_text, NULL);
    glCompileShader(fs);

    program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    mvp_location = glGetUniformLocation(program, "MVP");
    vpos_location = glGetAttribLocation(program, "vPos");
    vcol_location = glGetAttribLocation(program, "vCol");

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(vpos_location);
    glVertexAttribPointer(vpos_location, 2, GL_FLOAT, GL_FALSE,
        sizeof(Vertex), (void*)offsetof(Vertex, pos));

    glEnableVertexAttribArray(vcol_location);
    glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE,
        sizeof(Vertex), (void*)offsetof(Vertex, col));

    glEnable(GL_PROGRAM_POINT_SIZE);
}

void Renderer::tick(const Frame &frame) {

    if (!glfwWindowShouldClose(window)) {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        float ratio = width / (float)height;

        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT);

        glm::mat4 proj = glm::ortho(
            -ratio * zoom + cam_pos.x,
            ratio * zoom + cam_pos.x,
            -1.f * zoom + cam_pos.y,
            1.f * zoom + cam_pos.y,
            -1.f, 1.f
        );

        glUseProgram(program);
        glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(proj));

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        std::vector<Vertex> buffer;

        struct DrawCmd {
            GLenum mode;
            int start;
            int count;
        };

        std::vector<DrawCmd> commands;

        for (auto& pl : frame.planets) {
            int start, count;

            appendCircle(buffer,
                (float)pl.pos.x,
                (float)pl.pos.y,
                (float)pl.r,
                glm::vec3(0.2f, 0.6f, 1.f),
                start, count);

            commands.push_back({GL_TRIANGLE_FAN, start, count});
        }

        int p_start, p_count;
        appendPoints(buffer, frame.particles, p_start, p_count);

        commands.push_back({GL_POINTS, p_start, p_count});

        glBufferData(GL_ARRAY_BUFFER,
            buffer.size() * sizeof(Vertex),
            buffer.data(),
            GL_DYNAMIC_DRAW);

        for (auto& cmd : commands) {
            if (cmd.mode == GL_POINTS)
                glPointSize(frame.particle_r / zoom);

            glDrawArrays(cmd.mode, cmd.start, cmd.count);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    } else {
        glfwTerminate();
    }
}

void Renderer::quit() {

}