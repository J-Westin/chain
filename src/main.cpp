#include <iostream>
#include <thread>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/SVD>

#ifndef M_PI
    #define M_PI 3.14159265359
#endif // M_PI

const unsigned int N = 20;

const double     L = .1;
const double     g = 9.81;
const double gamma = g/L;
const double drawscale = 1./double(N);

typedef Eigen::Array<double, N, 1> ArrayNd;
typedef Eigen::Array<double, N, 2> xyblock;

struct jvec2{
public:
    double x, y;
    double length(){ return std::sqrt(x*x+y*y); }
    jvec2 operator-(){ return jvec2(-x, -y); }
};

#include <cmath>

double length(jvec2 u){
    return std::sqrt(u.x*u.x + u.y*u.y);
}

jvec2 rope_force()

xyblock posdifs(xyblock positions){
    xyblock difs = xyblock::Zero();

    difs(0,0) = positions(0,0);
    difs(0,1) = positions(0,1);

    for (unsigned int i = 1; i < N; i++){
        difs(i, 0) = positions(i, 0) - positions(i-1, 0);
        difs(i, 1) = positions(i, 1) - positions(i-1, 0);
    }

    return difs;
}

ArrayNd scalardifs(xyblock positions){
    xyblock difs = posdifs(positions);
    ArrayNd sdifs = ArrayNd::Zero();

}

//#include "jw/simple_chain.h"

#define GLEW_STATIC

#include <GL/glew.h>

#include <GLFW/glfw3.h>

double vertex_positions[2*(N+1)];


void fill_xy(double* position_array, ArrayNd thetas){
    double x_k(0.), y_k(0.);
    position_array[0] = x_k;
    position_array[1] = y_k;

    ArrayNd dx = drawscale*thetas.sin();
    ArrayNd dy = drawscale*thetas.cos();

    for (unsigned int i = 1; i <= N; i++){
        x_k += dx[i-1];
        y_k += dy[i-1];
        position_array[2*i  ] = x_k;
        position_array[2*i+1] = -y_k;
    }
    return;
}



int main(){

    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(600, 600, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    GLenum glew_init_error = glewInit();

    if ( glew_init_error != GLEW_OK ) {
        std::cerr << "GLEW initialization failed. Error code: " << glew_init_error << std::endl;
    }

    std::cout << "Using GL   version " << glGetString(GL_VERSION) << std::endl;
    std::cout << "Using GLSL version " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    GLuint vertex_positions_buffer_id;
    glGenBuffers(1, &vertex_positions_buffer_id);

    glBindBuffer(GL_ARRAY_BUFFER, vertex_positions_buffer_id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2*sizeof(double), 0);


    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);

        euler_integrate(thetas, omegas, 1./60.);

        fill_xy(vertex_positions, thetas);

        glBufferData(GL_ARRAY_BUFFER, (N+1)*2*sizeof(double), &vertex_positions, GL_DYNAMIC_DRAW);

        glDrawArrays(GL_LINE_STRIP, 0, N+1);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();

        std::this_thread::sleep_for(std::chrono::duration<double>(1./60.));
    }

    glfwTerminate();
    return 0;
}

