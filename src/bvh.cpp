#include "bvh.h"
#include "bvh_interface.h"
#include "common.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "stack"
#include <iostream>


// TODO Standard Feature
// Hierarchy traversal routine; you must implement this method and implement it carefully!
//
// The default implementation uses precompiled methods and only performs rudimentary updates to hitInfo.
// For correct normal interpolation, barycentric coordinates, etc., you have to implement this traversal yourself.
// If you have implemented the traversal method for assignment 4B, you might be able to reuse most of it.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
// - state;    current render state (containing scene, features, ...)
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
bool BVH::intersect(RenderState& state, Ray& ray, HitInfo& hitInfo) const
{
    if (!state.features.enableAccelStructure) {
        return intersectRayWithBVH(state, *this, ray, hitInfo);
    } else {
        // go through all spheres before going through meshes
        std::vector<Sphere> spheres = state.scene.spheres;
        for (int i = 0; i < spheres.size(); i++) {
            if (intersectRayWithShape(spheres[i], ray, hitInfo)) {
                hitInfo.material = spheres[i].material;
                return true;
            }
        }

        // create pointer to bvh and to the root node, initialize a stack with node pointers and push the root node, create a span of all the primitives.
        const BVHInterface* bvh = this;
        const BVHInterface::Node* currentNode = &bvh->nodes()[0];
        std::stack<const BVHInterface::Node*> stack;
        stack.push(currentNode);
        std::span<const BVHInterface::Primitive> primitives = bvh->primitives();

        // initialize intersected false, make 3 final intersection vectors, initialize the color to 0.1f(will be more and more green each recursive step)
        // intialize the final intersected node to null.
        bool intersected = false;
        Vertex final1;
        Vertex final2;
        Vertex final3;
        float color = 0.1f;
        const BVHInterface::Node* nodeKeep = nullptr;

        while (!stack.empty()) {
            // get the first node from the stack give it a wireframe and add color so next one would be more green.
            const BVHInterface::Node* node = stack.top();
            stack.pop();
            if (state.features.enableDebugDraw) {
                drawAABB(node->aabb, DrawMode::Wireframe, glm::vec3(1.0f - color, color, 0.0f), 1.0f);
                color += 0.01f;
                if (color > 1.0f) {
                    color = 1.0f;
                }
            }

            // if node is not leaf get the left and right children and add them to the stack so later can work on them.
            if (!node->isLeaf()) {
                const BVHInterface::Node* leftChild = &bvh->nodes()[node->leftChild()];
                const BVHInterface::Node* rightChild = &bvh->nodes()[node->rightChild()];
                float t = ray.t;

                auto leftBox = leftChild->aabb;
                auto rightBox = rightChild->aabb;
                // check if in aabb of left or right child (extra check)

                bool LeftIn = (ray.origin.x > leftBox.lower.x && ray.origin.x < leftBox.upper.x
                    && ray.origin.y > leftBox.lower.y && ray.origin.y < leftBox.upper.y
                    && ray.origin.z > leftBox.lower.z && ray.origin.z < leftBox.upper.z);
                bool RightIn = (ray.origin.x > rightBox.lower.x && ray.origin.x < rightBox.upper.x
                    && ray.origin.y > rightBox.lower.y && ray.origin.y < rightBox.upper.y
                    && ray.origin.z > rightBox.lower.z && ray.origin.z < rightBox.upper.z);
                // check if right or left to add to stack
                if (intersectRayWithShape(leftBox, ray) || LeftIn) {
                    stack.push(leftChild);
                    ray.t = t;
                }
                if (intersectRayWithShape(rightBox, ray) || RightIn) {
                    stack.push(rightChild);
                    ray.t = t;
                }
            }
            // if node is leaf iterate through its primitives and check if they intersect with triangle, if so update intersected,
            // nodeKeep, final vectors and hitInfo.
            else {
                uint32_t count = node->primitiveCount();
                uint32_t offset = node->primitiveOffset();
                for (uint32_t i = offset; i < offset + count; i++) {
                    BVHInterface::Primitive primitive = primitives[i];
                    if (intersectRayWithTriangle(primitive.v0.position, primitive.v1.position, primitive.v2.position, ray, hitInfo)) {
                        intersected = true;
                        // update material and barycentric coords
                        hitInfo.material = state.scene.meshes[primitive.meshID].material;
                        hitInfo.barycentricCoord = computeBarycentricCoord(primitive.v0.position, primitive.v1.position, primitive.v2.position, ray.origin + ray.direction * ray.t);
                        // update normal, use interpolated if it is on
                        if (state.features.enableNormalInterp) {
                            hitInfo.normal = interpolateNormal(primitive.v0.normal, primitive.v1.normal, primitive.v2.normal, hitInfo.barycentricCoord);
                        } else {
                            hitInfo.normal = glm::normalize(glm::cross(primitive.v1.position - primitive.v0.position, primitive.v2.position - primitive.v0.position));
                        }
                        // switch direction of normal if oposite from ray
                        if (glm::dot(ray.direction, hitInfo.normal) > 0.0f) {
                            hitInfo.normal = -hitInfo.normal;
                        }
                        if (state.features.enableDebugDraw) {
                            drawSphere((ray.direction * ray.t + ray.origin), 0.01f, glm::vec3(0.0f, 0.0f, 1.0f));
                        }
                        // update texture coords
                        hitInfo.texCoord = interpolateTexCoord(primitive.v0.texCoord, primitive.v1.texCoord, primitive.v2.texCoord, hitInfo.barycentricCoord);
                        // keep track of final intersection node
                        nodeKeep = node;
                        final1 = primitive.v0;
                        final2 = primitive.v1;
                        final3 = primitive.v2;
                    }
                }
            }
        }
        // if intersected draw the triangle, it intersected, make the wireframe of its aabb blue and put a little sphere at the point of intersection.
        if (intersected && state.features.enableDebugDraw) {
            if (nodeKeep != nullptr) {
                drawAABB(nodeKeep->aabb, DrawMode::Wireframe, glm::vec3(0.0f, 0.0f, 1.0f), 1.0f);
            }
            drawTriangle(final1, final2, final3);
            drawSphere((ray.direction * ray.t + ray.origin), 0.01f, glm::vec3(0.0f, 0.0f, 1.0f));
        }

        return intersected;
    }
}
