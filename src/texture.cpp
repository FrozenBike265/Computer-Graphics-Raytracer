#include "texture.h"
#include "render.h"
#include <framework/image.h>
#include <fmt/core.h>
#include <glm/gtx/string_cast.hpp>
#include <GLFW/glfw3.h>
#include <algorithm>
#include <iostream>

//extern RenderState gCurrentRenderState;

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    //return image.get_pixel(0);
    // we have to implement nearest texture sampling
    // we want to calculate the nearest pixel indices based on the texture coordinate
    // declarations
    int pxlOx, pxlOy, idx;

    // intializations
    pxlOx = static_cast<int>(texCoord.x * image.width);
    pxlOy = static_cast<int>(texCoord.y * image.height);

    // we have to clamp the pixel's coordinates to be sure
    // that we are inside the image boundaries
    pxlOx = std::max(0, std::min(pxlOx, image.width - 1));
    pxlOy = std::max(0, std::min(pxlOy, image.height - 1));

    // now we have to convert pixel coordinates to a linear index
    idx = pxlOy * image.width + pxlOx;
    /*
    if (gCurrentRenderState.features.enableTextureMapping && gCurrentRenderState.features.enableDebugDraw)
        displayDebugWindow(image, x, y);
    */
    // return the color of the pixel at idx - function already declared in image.h
    return image.get_pixel<3>(idx);
    
}
/*
void displayDebugWindow(const Image& Image, int sampledX, int sampledY)
{

}
*/



// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    //return image.get_pixel(0);

    // declarations
    int x, y, clidx1, clidx2, clidy1, clidy2;
    float fpx, fpy, dx, dy;
    glm::vec3 p00, p10, p01, p11, bil1, bil2, color;

    // first, we have to calculate the floating point pixel coordinates
    fpx = texCoord.x * image.width - 0.5f;
    fpy = texCoord.y * image.height - 0.5f;
    // now cast from float to int
    x = static_cast<int>(fpx);
    y = static_cast<int>(fpy);
    // find values between 0 and 1
    dx = fpx - floor(fpx);
    dy = fpy - floor(fpy);

    // we also have to make sure that the indices are not outside of bounds
    // to solve this we clamp the indices to be within the good range of the image
    clidx1 = std::clamp(x, 0, image.width - 1);
    clidx2 = std::clamp(x + 1, 0, image.width - 1);
    clidy1 = std::clamp(y, 0, image.height - 1);
    clidy2 = std::clamp(y + 1, 0, image.height - 1);

    // now we have get the four pixels that surround the image
    p00 = image.get_pixel<3>(clidy1 * image.width + clidx1);
    p10 = image.get_pixel<3>(clidy1 * image.width + clidx2);
    p01 = image.get_pixel<3>(clidy2 * image.width + clidx1);
    p11 = image.get_pixel<3>(clidy2 * image.width + clidx2);

    // now for the bilinear interpolation part we mix the pixels from the previous stepp
    bil1 = glm::mix(p00, p10, dx);
    bil2 = glm::mix(p01, p11, dx);
    color = glm::mix(bil1, bil2, dy);

    return color;
}