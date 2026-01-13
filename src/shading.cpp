#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>
#include <ranges>
#include <algorithm>
#include <GLFW/glfw3.h>
#include<vector>
#include<glm/glm.hpp>

// this function displays the visual debug for textures that are applied to the models
// the function shows a view of the texture space and how the geometry is mapped into it
// it also highlights with blue the ray intersection position
void displayTextureMapping(const Image& texture, const glm::vec2& texCoord)
{
    int idx, x, y;
    float dist;
    std::vector<unsigned char> pixelData(texture.width * texture.height * 3);
    glm::vec3 color;

    // i tried various ways to extract the pixeldata but the only one that I could make it work was this
    // implementation, which is similar to what we did in week 1
    for (y = 0; y < texture.height; y++)
        for (x = 0; x < texture.width; x++) {
            // get the color of each pixel
            color = texture.get_pixel<3>(y * texture.width + x);
            // get the index of each pixel
            idx = (y * texture.width + x) * 3;
            // we have to populate idx, idx+1, idx+2 because we work with RGB -> 3 pixels
            // get pixel is a predefined function from image.h which returns a value
            // from [0,1), this is why I multiply with 255
            pixelData[idx] = static_cast<unsigned char>(color.r * 255);
            pixelData[idx + 1] = static_cast<unsigned char>(color.g * 255);
            pixelData[idx + 2] = static_cast<unsigned char>(color.b * 255);
        }

    /*
    for (y = 0; y < texture.height; y++)
        for (x = 0; x < texture.width; x++)
        {
            color = texture.get_pixel<3>(y * texture.width + x);
            i = (y * texture.width + x) * 3;
            if (x / 10 == y / 10)
            {
                pixelData[i] = 255;
                pixelData[i + 1] = 255;
                pixelData[i + 2] = 255;
            }
            else
            {
                pixelData[i] = 0;
                pixelData[i + 1] = 0;
                pixelData[i + 2] = 0;
            }
        }
    */
    /*
    GLFWwindow* previous = glfwGetCurrentContext();

    GLFWwindow* window;
    if (!glfwInit())
        return;

    window = glfwCreateWindow(800, 800, "Texture Mapping Debug", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return;
    }

    glfwMakeContextCurrent(window);
    */

    // declare and initialize texture
    GLuint texId;
    glGenTextures(1, &texId);
    glBindTexture(GL_TEXTURE_2D, texId);

    // here we create an openGl texture with the image data
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture.width, texture.height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixelData.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // glBindTexture(GL_TEXTURE_2D, 0); //!!!!!
    glEnable(GL_TEXTURE_2D);

    // glBindTexture(GL_TEXTURE_2D, texId);//!!!!!

    // this is far enough to see very well the texture behind the teapot
    dist = 5;

    // we start rendering a quad with the texture applied to it
    glBegin(GL_QUADS);

    glTexCoord2f(0, 0);
    glVertex3f(-1, -1, dist);
    glTexCoord2f(1, 0);
    glVertex3f(1, -1, dist);
    glTexCoord2f(1, 1);
    glVertex3f(1, 1, dist);
    glTexCoord2f(0, 1);
    glVertex3f(-1, 1, dist);

    glEnd();
    // glBindTexture(GL_TEXTURE_2D, 0); //!!!!!
    glDisable(GL_TEXTURE_2D);

    // glDeleteTextures(1, &texId);

    // We highlight with the blue the texture coordinate
    glPointSize(10);
    glColor3f(0.0f, 0.0f, 1.0f);
    glBegin(GL_POINTS);
    glVertex3f(texCoord.x * 2 - 1, texCoord.y * 2 - 1, dist);
    glEnd();
}

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    glm::vec3 kdColor;

    // if texture, debugdraw and bvh of course are enabled, go to display the texture
    if (state.features.enableTextureMapping && state.features.enableDebugDraw && hitInfo.material.kdTexture) {
        displayTextureMapping(*hitInfo.material.kdTexture, hitInfo.texCoord);
    }

    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
        //if (state.features.enableDebugDraw)
            //displayTextureMappingDebug(*hitInfo.material.kdTexture, hitInfo.texCoord, kdColor);
    } else {
        //return hitInfo.material.kd;
        kdColor = hitInfo.material.kd;
    }
    return kdColor;
}



// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded color gradient. Feel free to modify this
    // static LinearGradient gradient = {
    //     .components = {
    //         { 0.1f, glm::vec3(215.f / 256.f, 210.f / 256.f, 203.f / 256.f) },
    //         { 0.22f, glm::vec3(250.f / 256.f, 250.f / 256.f, 240.f / 256.f) },
    //         { 0.5f, glm::vec3(145.f / 256.f, 170.f / 256.f, 175.f / 256.f) },
    //         { 0.78f, glm::vec3(255.f / 256.f, 250.f / 256.f, 205.f / 256.f) },
    //         { 0.9f, glm::vec3(170.f / 256.f, 170.f / 256.f, 170.f / 256.f) },
    //     }
    // };
    // TODO: come up with a good gradient
    static LinearGradient gradient = {
        .components = {
            { -1.0f, glm::vec3(1.0f, 0.1f, 0.1f) },
            // { 0.0f, glm::vec3(0.1f, 0.1f, 0.1f) },
            { 1.0f, glm::vec3(0.1f, 1.0f, 0.1f) }
        }
    };
    // static LinearGradient gradient = {
    //     .components = {
    //         { 0.0f, glm::vec3(1.0f, 0.1f, 0.1f) },
    //         { 0.5f, glm::vec3(0.1f, 0.1f, 1.0f) },
    //         { 1.0f, glm::vec3(0.1f, 1.0f, 0.1f) }
    //     }
    // };

    if (state.features.enableShading) {
        const glm::vec3 kd = sampleMaterialKd(state, hitInfo);
        switch (state.features.shadingModel) {
            case ShadingModel::Lambertian:
                return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
            case ShadingModel::Phong:
                return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
            case ShadingModel::BlinnPhong:
                return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
            case ShadingModel::LinearGradient:
                return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
            case ShadingModel::LinearGradientComparison:
                return computeLinearGradientModelComparison(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
            };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    std::vector<Component> sortedComponents;
    std::copy(components.cbegin(), components.cend(), std::back_inserter(sortedComponents));
    std::sort(sortedComponents.begin(), sortedComponents.end(), [](auto const& ca, auto const& cb) { return ca.t < cb.t; });
    int i = 0;
    for(; i < components.size(); i++) {
        if(ti <= components[i].t) {
            break;
        }
    }

    if(i == 0) {
        return components.front().color;
    }

    if(i == components.size()) {
        return components.back().color;
    }

    Component const& a = components[i - 1];
    Component const& b = components[i];
    float actual_t = (ti - a.t) / (b.t - a.t);

    return a.color * (1 - actual_t) + b.color * actual_t;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);
    return cos_theta * lightColor * gradient.sample(cos_theta);
}

glm::vec3 computeLinearGradientModelComparison(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);
    auto kd = sampleMaterialKd(state, hitInfo);
    auto phong = computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
    auto blinnPhong = computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
    auto d = glm::length(blinnPhong - phong);

    if(glm::dot(phong, phong) > glm::dot(blinnPhong, blinnPhong)) {
        return gradient.sample(d);
    } else {
        return gradient.sample(-d);
    }
}
