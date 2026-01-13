#include "extra.h"
#include "bvh.h"
#include "common.h"
#include "draw.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "texture.h"
#include <framework/trackball.h>
#include <GLFW/glfw3.h>

std::vector<glm::vec3> fp, lens, image, getfp;
std::vector<Ray> rays;

// I use this to handle the focus points
std::vector<glm::vec3> take_fp()
{
    return fp;
}

// I use this to handle lens
std::vector<glm::vec3> take_lens()
{
    return lens;
}

// I use this to handle image points
std::vector<glm::vec3> take_image()
{
    return image;
}

// I use this to handle the rays
std::vector<Ray> take_rays()
{
    return rays;
}

#include <algorithm>
#include <numeric>
#include <GL/gl.h>
#include <glm/gtx/polar_coordinates.hpp>

#include <functional>
#include <optional>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
    // 
    // Declarations
    int samples, x, y, sample;
    float paralel, lr, fd, fl, pi, k, planeDev, r, theta;
    long skipPoint, skipLens, pf, lf, rf;
    auto ef = features.extra;
    Ray ray, ray2;
    glm::vec2 pixelT, normPixelT;
    glm::vec3 cam, d, imagep, redPoint, lensPoint;
    glm::vec3 fpCenter, fpNormal, focalPoint, sColor;

    // Initializations
    skipPoint = 0; // Counter to determine when to debug focal points
    skipLens = 0; // Counter to determine when to debug lens points
    pf = 10; 
    lf = 15; // Frequency for adding focal points to the debug list
    rf = 100; // Frequency for adding camera rays to the debug list
    lr = ef.lr; // Lens radius
    fd = ef.fd; // Focal distance
    fl = ef.fl;  // Focal length
    samples = 8; // Number of rays sampled per pixel
    pi = 3.14159; // Pi value for trigonometric calculations
    cam = camera.position(); // Camera position
    d = camera.forward(); // Camera direction (view vector)
    fpCenter = cam + fd * d; // Focal plane center
    fpNormal = -d; // Focal plane normal (opposite to view direction)

    // Start of code
    for (y = 0;y<screen.resolution().y;y++)
        for (x = 0; x < screen.resolution().x; x++)
        {
            // Calculate the pixel position and normalize it to the[-1, 1] range
            pixelT = glm::vec2(x + 0.5f, y + 0.5f);
            normPixelT = glm::vec2(2.0f * pixelT.x / screen.resolution().x - 1.0f, 2.0f * pixelT.y / screen.resolution().y - 1.0f);
            
            // Generate a primary ray through the pixel
            ray = camera.generateRay(normPixelT);
            if (x % rf == 0 && y % rf == 0)
                rays.push_back(ray);

            // Calculate an intermediate point on the image plane
            imagep = ray.origin + 0.25f * ray.direction;
            image.push_back(imagep); // Save it for debugging purposes

            // Check if the ray intersects the focal plane
            paralel = glm::dot(ray.direction, fpNormal); // Dot product with the focal plane normal
            if (abs(paralel) < 1e-6f)
                continue; // Skip rays that are parallel to the focal plane

            // Calculate the intersection point on the focal plane
            k = glm::dot(fpCenter - ray.origin, fpNormal) / paralel;
            focalPoint = ray.origin + k * ray.direction;

            // Now let's start the debug part
            // Debug: Add focal points for visualization
            redPoint = glm::vec3(1.0f, 0.0f, 0.0f); // Red color for focal points
            skipPoint++;
            if (skipPoint % pf == 0)
                fp.push_back(focalPoint);

            // Check if the calculated focal point lies on the focal plane
            planeDev = glm::dot(focalPoint - fpCenter, fpNormal);
            if (abs(planeDev) > 0.001f)
                std::cout << "Point is not on plane"<<std::endl;
            
            // Initialize the pixel color for the current pixel
            glm::vec3 PixelColor(0.0f);

            // Generate secondary rays (lens rays) for depth of field effect
            for (sample = 0; sample < samples; sample++)
            {
                // Randomly sample a point on the lens (disk sampling)
                r = std::sqrt(static_cast<float>(rand()) / RAND_MAX) * lr * (fl / fd);
                theta = 2.0f * pi * static_cast<float>(rand()) / RAND_MAX;
                
                lensPoint = cam + (r * std::cos(theta) * camera.left()) + (r * std::sin(theta) * camera.up()) + (fl * d);
                
                // Debug: Add lens points for visualization
                skipLens++;
                if (skipLens % lf == 0)
                    lens.push_back(lensPoint);

                // Create a new ray from the lens point through the focal point
                ray2.origin = lensPoint;
                ray2.direction = glm::normalize(focalPoint - lensPoint);
                ray2.t = std::numeric_limits<float>::max(); // Set maximum ray distance

                // Render the ray and accumulate the color
                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(x + y * screen.resolution().x) }
                };
                sColor = renderRay(state, ray2, 0); // Recursive rendering
                PixelColor += sColor;
            }
            // Average the accumulated color and clamp it to [0, 1]
            PixelColor /= float(samples);
            screen.setPixel(x, y, glm::clamp(PixelColor, 0.0f, 1.0f));
        }
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

uint64_t binomialCoefficient(uint64_t n, uint64_t k) {
    uint64_t result = 1;
    for(int i = k + 1; i <= n; i++) {
        result *= i;
    }
    for(int i = 1; i <= n - k; i++) {
        result /= i;
    }
    return result;
}

void gaussianFilter(std::vector<glm::vec3>& pixels, Screen& screen, uint64_t filterSize) {
    std::vector<uint64_t> filter;
    std::vector<float> normalizedFilter;

    auto resolution = screen.resolution();

    uint64_t filterValuesSum = 0ull;
    for(int i = 0; i < 2 * filterSize + 1; i++) {
        filter.push_back(binomialCoefficient(2 * filterSize + 1 - 1, i));
        filterValuesSum += filter.back();
    }

    std::transform(filter.cbegin(), filter.cend(), std::back_inserter(normalizedFilter), [filterValuesSum](uint64_t a) {
        return static_cast<float>(a) / filterValuesSum;
    });

    std::vector<glm::vec3> auxBuffer(pixels.size(), glm::vec3(0.0f, 0.0f, 0.0f));

    // 1st pass - horizontal
    for(int y = 0; y < resolution.y; y++) {
        for(int x = 0; x < resolution.x; x++) {
            auto targetIndex = screen.indexAt(x, y);
            for(int i = 0; i < normalizedFilter.size(); i++) {
                int dx = i - filterSize;
                if(x + dx >= 0 && x + dx < resolution.x) {
                    auxBuffer[targetIndex] += normalizedFilter[i] * pixels[screen.indexAt(x + dx, y)];
                }
            }
        }
    }

    // 2nd pass - vertical
    for(int x = 0; x < resolution.x; x++) {
        for(int y = 0; y < resolution.y; y++) {
            auto targetIndex = screen.indexAt(x, y);
            for(int i = 0; i < normalizedFilter.size(); i++) {
                int dy = i - filterSize;
                if(y + dy >= 0 && y + dy < resolution.y) {
                    pixels[targetIndex] += normalizedFilter[i] * auxBuffer[screen.indexAt(x, y + dy)];
                }
            }
        }
    }
}

// use the formula used to calculate the relative luminance as used in sRGB
float luminance(glm::vec3 rgb) {
    return glm::dot(rgb, glm::vec3(0.2126f, 0.7152f, 0.0722f));
}

void grayScaleFilter(std::vector<glm::vec3>& pixels, Image& grayScaleImage, Screen& screen) {
    auto resolution = screen.resolution();
    auto& screenPixels = screen.pixels();
    auto filterWidth = grayScaleImage.width;
    auto filterHeight = grayScaleImage.height;
    std::vector<float> filter{};
    std::vector<glm::vec3> filteredPixels(pixels.size());

    for(int y = 0; y < filterHeight; y++) {
        for(int x = 0; x < filterWidth; x++) {
            filter.push_back(luminance(grayScaleImage.get_pixel(y * filterWidth + x)));
        }
    }

    float totalWeight = std::accumulate(filter.cbegin(), filter.cend(), 0.0f);
    std::transform(filter.cbegin(), filter.cend(), filter.begin(), [totalWeight](float f) { return f / totalWeight; });

    for (int y = 0; y < resolution.y; y++) {
        for (int x = 0; x < resolution.x; x++) {
            int targetIndex = screen.indexAt(x, y);

            for(int fy = 0; fy < filterHeight; fy++) {
                int dy = fy - filterHeight / 2;
                if(!(y + dy >= 0 && y + dy < resolution.y)) {
                    continue;
                }
                for(int fx = 0; fx < filterWidth; fx++) {
                    int dx = fx - filterWidth / 2;
                    if(!(x + dx >= 0 && x + dx < resolution.x)) {
                        continue;
                    }

                    filteredPixels[targetIndex] += pixels[screen.indexAt(x + dx, y + dy)] * filter[fy * filterWidth + fx];
                }
            }
        }
    }

    std::copy(filteredPixels.cbegin(), filteredPixels.cend(), pixels.begin());
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    float thresholdMin = features.extra.thresholdMin;
    float thresholdMax = features.extra.thresholdMax;
    BloomFunctionType bloomFunctionType = features.extra.bloomFunctionType;
    uint64_t filterSize = features.extra.bloomFilterSize;
    auto& pixels = image.pixels();

    std::vector<glm::vec3> thresholdedPixels {};
    thresholdedPixels.reserve(pixels.size());
    std::transform(pixels.cbegin(), pixels.cend(), std::back_inserter(thresholdedPixels), [thresholdMin, thresholdMax, bloomFunctionType](glm::vec3 const& rgb) -> glm::vec3 {
        float y = luminance(rgb);
        if (y < thresholdMin) {
            return glm::vec3(0.0f);
        }
        if (y >= thresholdMax) {
            return rgb;
        }
        if(thresholdMax == thresholdMin) {
            return rgb;
        }
        float x = (y - thresholdMin) / (thresholdMax - thresholdMin);
        float mask;

        switch(bloomFunctionType) {
            case CONSTANT:
                mask = 1.0f;
                break;
            case LINEAR:
                mask = x;
                break;
            case SQRT:
                mask = glm::sqrt(x);
                break;
        }

        return rgb * mask;
    });

    if(features.extra.bloomRenderType == MASK) {
        std::copy(thresholdedPixels.cbegin(), thresholdedPixels.cend(), pixels.begin());
        return;
    }


    if (features.extra.enableGrayScale && features.extra.bloomGrayscaleFilter) {
        grayScaleFilter(thresholdedPixels, *features.extra.bloomGrayscaleFilter, image);
    } else {
        gaussianFilter(thresholdedPixels, image, filterSize);
    }

    if(features.extra.bloomRenderType == CONTRIBUTION) {
        std::copy(thresholdedPixels.cbegin(), thresholdedPixels.cend(), pixels.begin());
        return;
    }

    for(int i = 0; i < pixels.size(); i++) {
        pixels[i] += thresholdedPixels[i];
    }
}

// Visual debug for:
// To verify the correctness of your implementation:
// Check whether the orthogonal basis created along the ray is correctly defined.
// Visualize the samples on the disk. Are the samples well (uniformly) distributed?
// Render the reflection ray bundle by drawing lines representing the secondary rays originating from a reflective surface point.
// Evaluate whether the bundle is centered around the perfect reflection direction.
// Check the spread of the bundle (i.e., narrower for higher shininess, wider for lower shininess) and explain your choice for how to reduce the size of your disk 

void allFourVerificationsVisualDebug(const glm::vec3& intersectPoint, const glm::vec3& perfectRefl, 
    const glm::vec3& basisU, const glm::vec3& basisV, const glm::vec3& basisW, float radius, int samples,
    const std::vector<glm::vec3>& sp, Sampler& sampler)
{
    //Declarations
    int circleSeg, i;
    float r, length, theta;
    glm::vec3 cp, ray;
    glm::vec3 sd, ed;

    // Initializations
    length = 0.5f;
    circleSeg = 64;

    // 1. Let's draw the orthogonal basis
    glBegin(GL_LINES);
    // Let's make U red
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(intersectPoint.x, intersectPoint.y, intersectPoint.z);
    glVertex3f(intersectPoint.x + length * basisU.x, intersectPoint.y + length * basisU.y, intersectPoint.z + length * basisU.z);
    // Let's make V green
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(intersectPoint.x, intersectPoint.y, intersectPoint.z);
    glVertex3f(intersectPoint.x + length * basisV.x, intersectPoint.y + length * basisV.y, intersectPoint.z + length * basisV.z);
    // Let's make W blue
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(intersectPoint.x, intersectPoint.y, intersectPoint.z);
    glVertex3f(intersectPoint.x + length * basisW.x, intersectPoint.y + length * basisW.y, intersectPoint.z + length * basisW.z);
    //end gl
    glEnd();

    // 2. Let's draw the samples on the disk (blue points)
    glColor3f(0.0f, 0.0f, 1.0f);
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    for (const auto& s : sp)
        glVertex3f(s.x, s.y, s.z);
    glEnd();

    // 3. Let's draw the sampling disk (green circle)
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_LINE_LOOP);
    for (i = 0; i < circleSeg; i++)
    {
        theta = 2.0f * glm::pi<float>() * float(i) / float(circleSeg);
        cp = intersectPoint + radius * (glm::cos(theta) * basisU + glm::sin(theta) * basisV);
        glVertex3f(cp.x, cp.y, cp.z);
    }
    glEnd();

    // 4. Draw the perfect reflection direction (pink line)
    glColor3f(1.0f, 0.0f, 1.0f);
    glBegin(GL_LINES);
    glVertex3f(intersectPoint.x, intersectPoint.y, intersectPoint.z);
    ray = intersectPoint + 2.0f * radius * perfectRefl;
    glVertex3f(ray.x, ray.y, ray.z);
    glEnd();

    // 5. Render ray bundle (while lines)
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    for (i = 0; i < samples; i++)
    {
        r = radius * sampler.next_1d();
        theta = 2.0f * glm::pi<float>() * sampler.next_1d();
        sd = glm::normalize(perfectRefl + r * (glm::cos(theta) * basisU + glm::sin(theta) * basisV));
        ed = intersectPoint + sd * 2.0f;
        glVertex3f(intersectPoint.x, intersectPoint.y, intersectPoint.z);
        glVertex3f(ed.x, ed.y, ed.z);
    }
    glEnd();
}

// uses the method from the book, finding the smallest component and changing it to 1.0f
glm::vec3 calculateBasis(glm::vec3 basisW)
{
    // finds smallest component in basisW
    float smallestComponent = std::abs(basisW[0]);
    int smallestComponentDirection = 0;
    if (smallestComponent > std::abs(basisW[1])) {
        smallestComponent = std::abs(basisW[1]);
        smallestComponentDirection = 1;
    }
    if (smallestComponent > std::abs(basisW[2])) {
        smallestComponent = std::abs(basisW[2]);
        smallestComponentDirection = 2;
    }

    // changes the smallest component to 1.0f
    glm::vec3 t;
    if (smallestComponentDirection == 0) {
        t = glm::vec3(1.0f, basisW.y, basisW.z);
    } else if (smallestComponentDirection == 1) {
        t = glm::vec3(basisW.x, 1.0f, basisW.z);
    } else {
        t = glm::vec3(basisW.x, basisW.y, 1.0f);
    }

    // if didnt work add 1 to basis
    if (glm::length(glm::cross(basisW, t)) < 1e-6f) {
        t += glm::vec3(1.0f, 0.0f, 0.0f);
    }

    // return perpendicular to this non colinear t
    return glm::normalize(glm::cross(basisW, t));
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...

    // we are going to use this for visual debug
    // we are going to store the sample points here
    std::vector<glm::vec3> sp;


    // if glossy not enabled
    if (!state.features.extra.enableGlossyReflection) {
        return;
    }

    // if not shiny or 0 ks return
    if (glm::length(hitInfo.material.ks) == 0.0f || hitInfo.material.shininess <= 0.0f) {
        return;
    }

    // get the reflected ray and the point of intersection
    Ray reflectedRay = generateReflectionRay(ray, hitInfo);
    glm::vec3 intersectPoint = ray.origin + ray.t * ray.direction;

    // change the radius based on the shinyness (inversly perportional)
    float radius = 1.0f / hitInfo.material.shininess;
    int samples = state.features.extra.numGlossySamples;
    glm::vec3 color = { 0.0f, 0.0f, 0.0f };

    // get basis perpendicular to reflected ray using fucntion calculateBasis
    glm::vec3 basisW = glm::normalize(reflectedRay.direction);
    glm::vec3 basisU = calculateBasis(basisW);
    glm::vec3 basisV = glm::normalize(glm::cross(basisW, basisU));

    //if (state.features.enableDebugDraw)
    //{
    //    debugBasis();
    //    debugDisk();
    //    debugRay();
    //}

    // render as many glossy rays as there are samples
    for (int i = 0; i < samples; i++) {
        // calculate a random distance from the origin of the disk
        float r = radius * state.sampler.next_1d();
        // calculate a random angle on the disk between 0 and 2 pi
        float theta = 2.0f * glm::pi<float>() * state.sampler.next_1d();

        // offset the origin of glossy to prevent self-intersections
        glm::vec3 originGlossy = intersectPoint + glm::normalize(reflectedRay.direction) * 0.001f;

        // direction of glossy is the reflected direction + the random offsets
        glm::vec3 directionGlossy = glm::normalize(reflectedRay.direction + r * glm::cos(theta) * basisU + r * glm::sin(theta) * basisV);
        Ray glossyRay = { originGlossy, directionGlossy };

        // when the ray is below the surface return a 0
        if (glm::dot(glossyRay.direction, hitInfo.normal) < 0.0f) {
            continue;
        }

        // add the sample point to use for visual debug
        sp.push_back(intersectPoint + r * (glm::cos(theta) * basisU + glm::sin(theta) * basisV));

        // add to color the result of renderRay of glossy
        color += renderRay(state, glossyRay, rayDepth + 1);
    }
    // add the avarage color of all the samples * ks to hitColor
    hitColor += (color / float(samples)) * hitInfo.material.ks;

    // use the renderstate and make an external method for all the visual debugs
    if (state.features.enableDebugDraw)
        allFourVerificationsVisualDebug(intersectPoint, reflectedRay.direction,
            basisU,basisV,basisW,radius,samples,sp,state.sampler);
}

std::pair<int, glm::vec2> getCubeMapCoordinates(glm::vec3 const& dir);
std::optional<glm::vec3> intersectRayWithAABBEnvMap(AxisAlignedBox const& aabb, Ray& ray);
using TextureSampler = glm::vec3 (*)(const Image&, const glm::vec2&);
glm::vec3 finiteCubeSampler(RenderState& state, Ray& ray, TextureSampler textureSampler);
glm::vec3 infiniteCubeSampler(RenderState& state, Ray& ray, TextureSampler textureSampler);

glm::vec3 infiniteSphereSampler(RenderState& state, Ray& ray, TextureSampler textureSampler);
void drawSquare(Image& texture, std::array<glm::vec3, 4> vertices);

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        float envMapSize = state.features.extra.cubeMapSize;

        // note: this gets *really* slow in ray tracing mode. Make sure to disable debug draw before rendering to a file!
        if(state.features.enableDebugDraw && !state.features.extra.useSphereMap) {
            drawSquare(*state.scene.cubeMap[0], std::array<glm::vec3, 4>{
                glm::vec3{ 1.0f, -1.0f, -1.0f } * envMapSize,
                glm::vec3{ 1.0f, -1.0f, 1.0f } * envMapSize,
                glm::vec3{ 1.0f, 1.0f, 1.0f } * envMapSize,
                glm::vec3{ 1.0f, 1.0f, -1.0f } * envMapSize
            });
            drawSquare(*state.scene.cubeMap[1], std::array<glm::vec3, 4>{
                glm::vec3{ -1.0f, -1.0f, 1.0f } * envMapSize,
                glm::vec3{ -1.0f, -1.0f, -1.0f } * envMapSize,
                glm::vec3{ -1.0f, 1.0f, -1.0f } * envMapSize,
                glm::vec3{ -1.0f, 1.0f, 1.0f } * envMapSize
            });
            drawSquare(*state.scene.cubeMap[2], std::array<glm::vec3, 4>{
                glm::vec3{ 1.0f, 1.0f, 1.0f } * envMapSize,
                glm::vec3{ -1.0f, 1.0f, 1.0f } * envMapSize,
                glm::vec3{ -1.0f, 1.0f, -1.0f } * envMapSize,
                glm::vec3{ 1.0f, 1.0f, -1.0f } * envMapSize
            });
            drawSquare(*state.scene.cubeMap[3], std::array<glm::vec3, 4>{
                glm::vec3{ 1.0f, -1.0f, -1.0f } * envMapSize,
                glm::vec3{ -1.0f, -1.0f, -1.0f } * envMapSize,
                glm::vec3{ -1.0f, -1.0f, 1.0f } * envMapSize,
                glm::vec3{ 1.0f, -1.0f, 1.0f } * envMapSize
            });
            drawSquare(*state.scene.cubeMap[4], std::array<glm::vec3, 4>{
                glm::vec3{ 1.0f, -1.0f, 1.0f } * envMapSize,
                glm::vec3{ -1.0f, -1.0f, 1.0f } * envMapSize,
                glm::vec3{ -1.0f, 1.0f, 1.0f } * envMapSize,
                glm::vec3{ 1.0f, 1.0f, 1.0f } * envMapSize
            });
            drawSquare(*state.scene.cubeMap[5], std::array<glm::vec3, 4>{
                glm::vec3{ -1.0f, -1.0f, -1.0f } * envMapSize,
                glm::vec3{ 1.0f, -1.0f, -1.0f } * envMapSize,
                glm::vec3{ 1.0f, 1.0f, -1.0f } * envMapSize,
                glm::vec3{ -1.0f, 1.0f, -1.0f } * envMapSize
            });
        }

        TextureSampler textureSampler;

        if(state.features.enableBilinearTextureFiltering) {
            textureSampler = sampleTextureBilinear;
        } else {
            textureSampler = sampleTextureNearest;
        }

        glm::vec3 color;

        if(state.features.extra.useSphereMap) {
            color = infiniteSphereSampler(state, ray, textureSampler);
        } else {
            if(glm::epsilonEqual(state.features.extra.cubeMapSize, 50.0f, 1e-4f)) {
                color = infiniteCubeSampler(state, ray, textureSampler);
            } else {
                color = finiteCubeSampler(state, ray, textureSampler);
            }
        }

        if(state.features.enableDebugDraw) {
            drawRay(ray, color);
        }

        return color;
    } else {
        return glm::vec3(0.f);
    }
}

glm::vec3 finiteCubeSampler(RenderState& state, Ray& ray, TextureSampler textureSampler) {
    float envMapSize = state.features.extra.cubeMapSize;
    AxisAlignedBox aabb = { { -envMapSize, -envMapSize, -envMapSize }, { envMapSize, envMapSize, envMapSize } };
    auto res = intersectRayWithAABBEnvMap(aabb, ray);

    if(!res) {
        return glm::vec3(0.0f);
    }

    auto pos = *res;
    constexpr float eps = 1e-4f;

    if(glm::epsilonEqual(pos.x, envMapSize, eps)) {
        auto u = (pos.z + envMapSize) / (2.0f * envMapSize);
        auto v = (pos.y + envMapSize) / (2.0f * envMapSize);
        return textureSampler(*state.scene.cubeMap[0], { u, 1.0f - v });
    } else if(glm::epsilonEqual(pos.x, -envMapSize, eps)) {
        auto u = 1.0f - (pos.z + envMapSize) / (2.0f * envMapSize);
        auto v = (pos.y + envMapSize) / (2.0f * envMapSize);
        return textureSampler(*state.scene.cubeMap[1], { u, 1.0f - v });
    } else if(glm::epsilonEqual(pos.y, envMapSize, eps)) {
        auto u = 1.0f - (pos.x + envMapSize) / (2.0f * envMapSize);
        auto v = 1.0f - (pos.z + envMapSize) / (2.0f * envMapSize);
        return textureSampler(*state.scene.cubeMap[2], { u, 1.0f - v });
    } else if(glm::epsilonEqual(pos.y, -envMapSize, eps)) {
        auto u = 1.0f - (pos.x + envMapSize) / (2.0f * envMapSize);
        auto v = (pos.z + envMapSize) / (2.0f * envMapSize);
        return textureSampler(*state.scene.cubeMap[3], { u, 1.0f - v });
    } else if(glm::epsilonEqual(pos.z, envMapSize, eps)) {
        auto u = 1.0f - (pos.x + envMapSize) / (2.0f * envMapSize);
        auto v = (pos.y + envMapSize) / (2.0f * envMapSize);
        return textureSampler(*state.scene.cubeMap[4], { u, 1.0f - v });
    } else if(glm::epsilonEqual(pos.z, -envMapSize, eps)) {
        auto u = (pos.x + envMapSize) / (2.0f * envMapSize);
        auto v = (pos.y + envMapSize) / (2.0f * envMapSize);
        return textureSampler(*state.scene.cubeMap[5], { u, 1.0f - v });
    } else {
        // SHOULD NEVER HAPPEN
        return glm::vec3(0.0f);
    }
}

glm::vec3 infiniteCubeSampler(RenderState& state, Ray& ray, TextureSampler textureSampler) {
    auto [idx, uv] = getCubeMapCoordinates(ray.direction);
    return textureSampler(*state.scene.cubeMap[idx], uv);
}

glm::vec3 infiniteSphereSampler(RenderState& state, Ray& ray, TextureSampler textureSampler) {
    auto p = glm::polar(ray.direction);
    auto u = p.y / (2 * glm::pi<float>()) + 0.5f;
    auto v = p.x / glm::pi<float>() + 0.5f;
    auto uv = glm::vec2(u, 1.0f - v);
    return textureSampler(*state.scene.sphereMap, uv);
}

void drawSquare(Image& texture, std::array<glm::vec3, 4> vertices) {
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture.width, texture.height, 0, GL_RGB, GL_UNSIGNED_BYTE, texture.get_data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 1.0f);
        glVertex3fv(glm::value_ptr(vertices[0]));

        glTexCoord2f(1.0f, 1.0f);
        glVertex3fv(glm::value_ptr(vertices[1]));

        glTexCoord2f(1.0f, 0.0f);
        glVertex3fv(glm::value_ptr(vertices[2]));

        glTexCoord2f(0.0f, 0.0f);
        glVertex3fv(glm::value_ptr(vertices[3]));
    glEnd();
    glDisable(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, 0);
    glDeleteTextures(1, &tex);
}

// based on description from chapter 20.4 of "Computer graphics : principles and practice", 3rd edition, ISBN 9780133373721
std::pair<int, glm::vec2> getCubeMapCoordinates(glm::vec3 const& dir) {
    auto x = dir.x;
    auto y = dir.y;
    auto z = dir.z;
    auto magX = glm::abs(x);
    auto magY = glm::abs(y);
    auto magZ = glm::abs(z);
    auto maxAxis = glm::max(glm::max(magX, magY), magZ);
    int face = 0;
    float u = 0;
    float v = 0;

    if(maxAxis == magX) {
        if(x <= 0) {
            face = 1;
            u = -z;
        } else {
            face = 0;
            u = z;
        }
        v = y;
    } else if(maxAxis == magY) {
        if(y <= 0) {
            face = 3;
            v = z;
        } else {
            face = 2;
            v = -z;
        }
        u = -x;
    } else if(maxAxis == magZ) {
        if(z <= 0) {
            face = 5;
            u = x;
        } else {
            face = 4;
            u = -x;
        }
        v = y;
    }

    auto t = 2 * maxAxis;
    return std::pair<int, glm::vec2>{ face, glm::vec2(u / t + 0.5f, 1.0f - (v / t + 0.5f)) };
}

std::optional<glm::vec3> intersectRayWithAABBEnvMap(AxisAlignedBox const& aabb, Ray& ray) {
    // I copied this code from my implementation of assignment 4b
    float txmin = (aabb.lower.x - ray.origin.x) / ray.direction.x;
    float txmax = (aabb.upper.x - ray.origin.x) / ray.direction.x;
    float tinx = glm::min(txmin, txmax);
    float toutx = glm::max(txmin, txmax);
    float tymin = (aabb.lower.y - ray.origin.y) / ray.direction.y;
    float tymax = (aabb.upper.y - ray.origin.y) / ray.direction.y;
    float tiny = glm::min(tymin, tymax);
    float touty = glm::max(tymin, tymax);
    float tzmin = (aabb.lower.z - ray.origin.z) / ray.direction.z;
    float tzmax = (aabb.upper.z - ray.origin.z) / ray.direction.z;
    float tinz = glm::min(tzmin, tzmax);
    float toutz = glm::max(tzmin, tzmax);

    float lastIn = glm::max(glm::max(tinx, tiny), tinz);
    float firstOut = glm::min(glm::min(toutx, touty), toutz);

    if(lastIn > firstOut) {
        // when ray is cast from outside of the cube map and doesn't ever hit it
        return std::nullopt;
    }

    if(lastIn > 0) {
        // when the ray is cast from outside of the cube map and hits it
        ray.t = lastIn;
    } else {
        // when the ray is cast from inside of the cube map
        ray.t = firstOut;
    }

    return { ray.origin + ray.direction * ray.t };
}