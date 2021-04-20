#ifndef SCENE_H
#define SCENE_H

/**
 * @brief Class to manage the scene for rendering
 */

#include "camera.h"
#include "meshrenderer.h"
#include "lights/light.h"
#include "datatypes/aabox.h"

#include <vector>

class Scene
{
public:
    Scene();

    Camera* getCam();
    std::shared_ptr<MeshRenderer> addPrimitive(const std::string& mesh_name, std::shared_ptr<Material>& mat);
    std::shared_ptr<MeshRenderer> addPrimitive(const std::shared_ptr<MeshRenderer>& mesh);
    void addLight(const std::shared_ptr<Light>& light);
    AABox getBoundingBox() const;

    const std::vector<std::shared_ptr<MeshRenderer>>& getPrimitives() const;
    const std::vector<std::shared_ptr<Light>>& getLights() const;

    static Scene* mainScene;

protected:
    Camera m_cam;
    std::vector<std::shared_ptr<MeshRenderer>> m_objects;
    std::vector<std::shared_ptr<Light>> m_lights;
};

#endif // SCENE_H
