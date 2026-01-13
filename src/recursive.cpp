#include "recursive.h"
#include "bvh_interface.h"
#include "draw.h"
#include "extra.h"
#include "intersect.h"
#include "light.h"

// This function is provided as-is. You do not have to implement it.
// Given a range of rays, render out all rays and average the result
glm::vec3 renderRays(RenderState& state, std::span<const Ray> rays, int rayDepth)
{
    glm::vec3 L { 0.f };
    for (const auto& ray : rays) {
        L += renderRay(state, ray, rayDepth);
    }
    return L / static_cast<float>(rays.size());
}

// This method is provided as-is. You do not have to implement it.
// Given a camera ray (or secondary ray), tests for a scene intersection, and
// dependent on the results, evaluates the following functions which you must
// implement yourself:
// - computeLightContribution() and its submethods
// - renderRaySpecularComponent(), renderRayTransparentComponent(), renderRayGlossyComponent()
glm::vec3 renderRay(RenderState& state, Ray ray, int rayDepth)
{
    // Trace the ray into the scene. If nothing was hit, return early
    HitInfo hitInfo;
    if (!state.bvh.intersect(state, ray, hitInfo)) {
        if (state.features.enableDebugDraw && !state.features.extra.enableEnvironmentMap) {
            drawRay(ray, glm::vec3(1, 0, 0));
        }
        return sampleEnvironmentMap(state, ray);
    }

    // Return value: the light along theray
    // Given an intersection, estimate the contribution of scene lights at this intersection
    glm::vec3 Lo = computeLightContribution(state, ray, hitInfo);

    // Draw an example debug ray for the incident ray (feel free to modify this for yourself)
    drawRay(ray, glm::vec3(1.0f));

    // Given that recursive components are enabled, and we have not exceeded maximum depth,
    // estimate the contribution along these components
    if (rayDepth < 6) {
        bool isReflective = glm::any(glm::notEqual(hitInfo.material.ks, glm::vec3(0.0f)));
        bool isTransparent = hitInfo.material.transparency != 1.f;

        // Default, specular reflections
        if (state.features.enableReflections && !state.features.extra.enableGlossyReflection && isReflective) {
            renderRaySpecularComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Alternative, glossy reflections
        if (state.features.enableReflections && state.features.extra.enableGlossyReflection && isReflective) {
            renderRayGlossyComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Transparency passthrough
        if (state.features.enableTransparency && isTransparent) {
            renderRayTransparentComponent(state, ray, hitInfo, Lo, rayDepth);
        }
    }

    return Lo;
}
// TODO: Standard feature
// Given an incident ray and a intersection point, generate a mirrored ray
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a reflected ray
// This method is unit-tested, so do not change the function signature.
Ray generateReflectionRay(Ray ray, HitInfo hitInfo)
{
    // TODO: generate a mirrored ray
    //       if you use glm::reflect, you will not get points for this method!
    glm::vec3 direction = glm::normalize(ray.direction);
    // Get direction of reflected ray
    glm::vec3 resultDirection = glm::normalize(direction - 2.0f * (glm::dot(direction, hitInfo.normal)) * hitInfo.normal);
    glm::vec3 origin = ray.origin + ray.t * ray.direction;
    // return new ray with slight offset for intersection
    Ray result = { origin + resultDirection * 0.001f, resultDirection };
    return result;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a passthrough ray for transparency,
// starting at the intersection point and continuing in the same direction.
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a passthrough ray for transparency
// This method is unit-tested, so do not change the function signature.
Ray generatePassthroughRay(Ray ray, HitInfo hitInfo)
{
    // TODO: generate a passthrough ray
    glm::vec3 origin = ray.origin + ray.t * ray.direction;
    // return new ray with slight offset for intersection
    return { origin + 0.001f * ray.direction, ray.direction };
}
// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a mirrored ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and adding the result times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRaySpecularComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generateReflectionRay()
    // get reflected ray
    Ray reflectedRay = generateReflectionRay(ray, hitInfo);
    // recursivley get reflected ray color
    glm::vec3 color = renderRay(state, reflectedRay, rayDepth + 1);
    // add the color * the ks
    hitColor += color * hitInfo.material.ks;

    // if we are in the recursive tree at the first level draw a blue ray (initial output ray, if second level first reflection is red
    // every level lower is slightly more green and blue mix and less red (pink)
    if (state.features.enableDebugDraw) {
        if (rayDepth == 0) {
            drawRay(ray, glm::vec3(0.0f, 0.0f, 1.0f));
        }
        if (rayDepth == 1) {
            drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        }
        if (rayDepth == 2) {
            drawRay(ray, glm::vec3(0.8f, 0.2f, 0.2f));
        }
        if (rayDepth == 3) {
            drawRay(ray, glm::vec3(0.6f, 0.4f, 0.4f));
        }
        if (rayDepth == 4) {
            drawRay(ray, glm::vec3(0.4f, 0.6f, 0.6f));
        }
        if (rayDepth == 5) {
            drawRay(ray, glm::vec3(0.2f, 0.8f, 0.8f));
            drawRay(reflectedRay, glm::vec3(0.0f, 1.0f, 1.0f));
        }
        glm::vec3 origin = ray.origin + ray.t * ray.direction;
        // draw normal for each
        drawRay(Ray { origin, hitInfo.normal, 0.05f }, glm::vec3(0.0f, 1.0f, 0.0f));
    }
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a passthrough transparent ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and correctly alpha blending the result with the current intersection's hit color
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRayTransparentComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generatePassthroughRay()
    float transparency = hitInfo.material.transparency;

    // get passThrough and call renderRay with it for color
    Ray passthroughRay = generatePassthroughRay(ray, hitInfo);
    glm::vec3 transmittedColor = renderRay(state, passthroughRay, rayDepth + 1);
    // set hit color *  alpha and color after hit * 1- alpha
    hitColor = (transparency)*hitColor + (1 - transparency) * transmittedColor;
    if (state.features.enableDebugDraw) {
        // draw ray and sphere of with hit color
        drawRay(ray, hitColor);
        glm::vec3 origin = ray.origin + ray.t * ray.direction;
        drawSphere(origin, 0.01f, hitColor);

        Ray normal = { origin, hitInfo.normal, 0.05f };

        Ray normalExtended = { normal.direction * normal.t + normal.origin, hitInfo.normal, 0.05f };
        // draw the first part of the normal ray as the kd value and the second as the transparency value
        drawRay(normal, hitInfo.material.kd);
        drawRay(normalExtended, { transparency, transparency, transparency });
    }
}