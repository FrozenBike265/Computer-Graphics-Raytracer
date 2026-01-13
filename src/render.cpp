#include "render.h"
#include "bvh_interface.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include "draw.h"
#include "extra.h"
#include "light.h"
#include "recursive.h"
#include "sampler.h"
#include "screen.h"
#include "shading.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

// This function is provided as-is. You do not have to implement it.
// Given relevant objects (scene, bvh, camera, etc) and an output screen, multithreaded fills
// each of the pixels using one of the below `renderPixel*()` functions, dependent on scene
// configuration. By default, `renderPixelNaive()` is called.
void renderImage(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    // Either directly render the image, or pass through to extra.h methods
    if (features.extra.enableDepthOfField) {
        renderImageWithDepthOfField(scene, bvh, features, camera, screen);
    } else if (features.extra.enableMotionBlur) {
        renderImageWithMotionBlur(scene, bvh, features, camera, screen);
    } else {
#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
        for (int y = 0; y < screen.resolution().y; y++) {
            for (int x = 0; x != screen.resolution().x; x++) {
                // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
                // Note; we seed the sampler for consistenct behavior across frames
                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
                };
                auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
                auto L = renderRays(state, rays);
                screen.setPixel(x, y, L);
            }
        }
    }

    // Pass through to extra.h for post processing
    if (features.extra.enableBloomEffect) {
        postprocessImageWithBloom(scene, features, camera, screen);
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples for this pixel.
// This method forwards to `generatePixelRaysMultisampled` and `generatePixelRaysStratified` when necessary.
std::vector<Ray> generatePixelRays(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    if (state.features.numPixelSamples > 1) {
        if (state.features.enableJitteredSampling) {
            return generatePixelRaysStratified(state, camera, pixel, screenResolution);
        } else {
            return generatePixelRaysUniform(state, camera, pixel, screenResolution);
        }
    } else {
        // Generate single camera ray placed at the pixel's center
        // Note: (-1, -1) at the bottom left of the screen,
        //       (+1, +1) at the top right of the screen.
        glm::vec2 position = (glm::vec2(pixel) + 0.5f) / glm::vec2(screenResolution) * 2.f - 1.f;
        return { camera.generateRay(position) };
    }
}

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// uniformly throughout this pixel.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is unit-tested, so do not change the function signature.
std::vector<Ray> generatePixelRaysUniform(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples camera rays uniformly distributed across the pixel. Use
    // Hint; use `state.sampler.next*d()` to generate random samples in [0, 1).
    // auto numSamples = state.features.numPixelSamples;
    // std::vector<Ray> rays;
    // ...
    // return rays;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //state.sampler.seed(seed);
    //Sampler state.sampler(seed);

    // my code has 2 parts:
    // 1. make the corresponding rays for the function
    // 2. visual debug for the rays

    // start of first part

    // declarations
    float step, u, v, pp;
    glm::vec2 sp, np;
    int i, numSamples = state.features.numPixelSamples;
    std::vector<Ray> rays;
    std::vector<glm::vec2> uvCoordinates;

    // i use this for sample distribution within the pixel
    step = 1.0f / std::sqrt(numSamples);
    
    // generate the uniform rays
    for (i = 0; i < numSamples; i++) 
    {
        // used the hint for random offset
        u = state.sampler.next_1d();
        // used the hint for random offset
        v = state.sampler.next_1d();
        sp = glm::vec2(pixel) + glm::vec2(u, v);
        np = (sp / glm::vec2(screenResolution)) * 2.0f - 1.0f;
        rays.push_back(camera.generateRay(np));
        // I use uvCoordinates to deal with the grid in visual debug
        uvCoordinates.push_back(glm::vec2(u, v));
    }
    // grid size
    float gridSize = int(std::sqrt(numSamples));

    // here part 2 starts

    // If debugDraw is enabled:
    if (state.features.enableDebugDraw)
    {
        // I use previous so that after you close the pop up window, you can go back to the main application
        GLFWwindow* previous = glfwGetCurrentContext();
        if (!glfwInit())
        {
            return rays;
        }
        GLFWwindow* window = glfwCreateWindow(500, 500, "Visual Debug Uniform Distribution", nullptr, nullptr);
        if (!window)
        {
            glfwTerminate();
            return rays;
        }
        // Make the current window the one to use
        glfwMakeContextCurrent(window);
        // I made the one from the BOOK
        // So, here we have white color
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glEnable(GL_POINT_SMOOTH);
        // Bigger points for better visualization
        glPointSize(10.0f);

        while (!glfwWindowShouldClose(window))
        {
            glClear(GL_COLOR_BUFFER_BIT);

            // Draw the grid
            glColor3f(0.8f, 0.8f, 0.8f);
            glBegin(GL_LINES);

            for (i = 0; i <= gridSize; i++) 
            {
                pp = i / gridSize;
                glVertex2f(pp - 0.5f, -0.5f);
                glVertex2f(pp - 0.5f, 0.5f);
                glVertex2f(-0.5f, pp - 0.5f);
                glVertex2f(0.5f, pp - 0.5f);
            }
            glEnd();

            // Draw the points for each ray
            glBegin(GL_POINTS);
            // Use blue color as in the BOOK
            glColor3f(0.0f, 0.0f, 1.0f);

            for (const auto& coord : uvCoordinates)
            {
                glVertex2f(coord.x - 0.5f, coord.y - 0.5f);
            }

            glEnd();
            glfwSwapBuffers(window);
            glfwPollEvents();
        }
        glfwDestroyWindow(window);
        glfwMakeContextCurrent(previous);
        //glfwTerminate();
    }
    
    return rays;
}

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// using jittered sampling throughout this pixel. Given NxN cells across one pixel, each ray sample is randomly
// placed somewhere within a cell.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
std::vector<Ray> generatePixelRaysStratified(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples * numSamples camera rays as jittered samples across the pixel.
    // Hint; use `state.sampler.next*d()` to generate random samples in [0, 1).
    //auto numSamples = static_cast<uint32_t>(std::round(std::sqrt(float(state.features.numPixelSamples))));
    //std::vector<Ray> rays;
    //...
    //return rays;

    // my code has 2 parts:
    // 1. make the corresponding rays for the function
    // 2. visual debug for the rays

    // start of first part

    // declarations

    float step, ox, oy, u, v, pp;
    glm::vec2 sp, np;
    std::vector<glm::vec2> uvCoordinates; //!!!!!!!!!!!!!!!
    int i, j, numSamples = static_cast<int>(std::sqrt(state.features.numPixelSamples));
    std::vector<Ray> rays;

    rays.reserve(numSamples * numSamples);

    // I use this for sample distribution within the pixel
    step = 1.0f / numSamples;

    // generate the jittered rays
    for (i = 0; i < numSamples; i++)
        for (j = 0; j < numSamples; j++) {
            // now, we generate random offsets within each subpixel region
            ox = state.sampler.next_1d() * step;
            oy = state.sampler.next_1d() * step;

            // calculate the positions horizontal and vertical
            u = (j + ox) * step;
            v = (i + oy) * step;

            sp = glm::vec2(pixel) + glm::vec2(u, v);
            np = (sp / glm::vec2(screenResolution)) * 2.0f - 1.0f;

            // generate the rays and store them
            rays.push_back(camera.generateRay(np));
            // store in uvCoordinates to use for visual debug
            uvCoordinates.push_back(glm::vec2(u, v));
        }

    // here part 2 starts
    // this is the same as in uniform distribution

    // If debugDraw is enabled:
    if (state.features.enableDebugDraw) {
        // I use previous so that after you close the pop up window, you can go back to the main application
        GLFWwindow* previous = glfwGetCurrentContext();
        if (!glfwInit()) {
            return rays;
        }
        GLFWwindow* window = glfwCreateWindow(500, 500, "Visual Debug Stratified Rays", nullptr, nullptr);
        if (!window) {
            glfwTerminate();
            return rays;
        }
        // Make the current window the one to use
        glfwMakeContextCurrent(window);
        // I made the one from the BOOK
        // So, here we have white color
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glEnable(GL_POINT_SMOOTH);
        // Bigger points for better visualization
        glPointSize(10.0f);

        while (!glfwWindowShouldClose(window)) {
            glClear(GL_COLOR_BUFFER_BIT);

            // Draw the grid
            glColor3f(0.8f, 0.8f, 0.8f);
            glBegin(GL_LINES);
            for (i = 0; i <= numSamples; i++) {
                pp = i / float(numSamples);
                glVertex2f(pp - 0.5f, -0.5f);
                glVertex2f(pp - 0.5f, 0.5f);
                glVertex2f(-0.5f, pp - 0.5f);
                glVertex2f(0.5f, pp - 0.5f);
            }
            glEnd();
            // Draw the points for each ray
            glBegin(GL_POINTS);
            glColor3f(0.0f, 0.0f, 1.0f);
            // Use blue color as in the BOOK
            for (const auto& coord : uvCoordinates) {
                glVertex2f(coord.x - 0.5f, coord.y - 0.5f);
            }

            glEnd();
            glfwSwapBuffers(window);
            glfwPollEvents();
        }
        glfwDestroyWindow(window);
        glfwMakeContextCurrent(previous);
        // glfwTerminate();
    }

    return rays;
}


