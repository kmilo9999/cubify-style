#include "scene.h"

Scene* Scene::main = nullptr;

Scene::Scene()
{
    if (main == nullptr) {
        main = this;
    }
}

Camera* Scene::getCam()
{
    return &m_cam;
}

const std::vector<std::shared_ptr<MeshRenderer>>& Scene::getPrimitives() const
{
    return m_objects;
}

const std::vector<std::shared_ptr<Light>>& Scene::getLights() const
{
    return m_lights;
}

std::shared_ptr<MeshRenderer> Scene::addPrimitive(const std::string &mesh_name, std::shared_ptr<Material>& mat)
{
    // first, try to create a mesh
    std::shared_ptr<MeshRenderer> primitive = std::make_shared<MeshRenderer>();
    if (!primitive->load_mesh(mesh_name)) {
        // failed, return nullptr
        return nullptr;
    }
    // else, add to the scene and return the pointer
    primitive->set_material(mat);
    m_objects.emplace_back(primitive);
    return primitive;
}

std::shared_ptr<MeshRenderer> Scene::addPrimitive(const std::shared_ptr<MeshRenderer> &mesh)
{
    m_objects.emplace_back(mesh);
    return mesh;
}

void Scene::addLight(const std::shared_ptr<Light> &light)
{
    m_lights.emplace_back(light);
}

AABox Scene::getBoundingBox() const
{
    AABox ret;
    for (const auto ptr : m_objects) {
        ret.makeContain(ptr->get_boundingbox());
    }
    return ret;
}
