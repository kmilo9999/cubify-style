#include "demoscene.h"
#include "graphics/materials/toonmaterial.h"

DemoScene::DemoScene() : Scene()
{
    // Scene for the demo
    // add sphere
    std::shared_ptr<ToonMaterial> toon = std::make_shared<ToonMaterial>();

    std::shared_ptr<ToonMaterial> toon_bunny = std::make_shared<ToonMaterial>();
    toon_bunny->diffuseColor = Eigen::Vector3f(0.8f, 0.8f, 0.3f);

    std::shared_ptr<Material> mat = std::dynamic_pointer_cast<Material>(toon);
    addPrimitive("./meshes/sphere.obj", mat);

    // mat = std::dynamic_pointer_cast<Material>(toon_bunny);
    // addPrimitive("./meshes/bunny.obj", mat)->set_location(Eigen::Vector3f(0, 2.f, 0));

    // add lighting
    std::shared_ptr<DirectionalLight> lit = std::make_shared<DirectionalLight>(Eigen::Vector3f(0.5f, 0.5f, 0.65f), Eigen::Vector3f(-1.f, -1.f, 0.f));
    addLight(lit);

    // add camera
    Camera* cam = getCam();
    cam->setPosition(Eigen::Vector3f(0, 0, 5));
    cam->lookAt(Eigen::Vector3f(0, 2, -5), Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(0, 1, 0));
    // cam->setTarget(Eigen::Vector3f(0, 2, 0));

    // perspect will set in the rendering pipeline.
}

void DemoScene::getCubifyMeshes(std::vector<std::shared_ptr<Mesh> > &meshes)
{
    meshes.clear();
    for (size_t i = 0; i < m_objects.size(); ++i) {
        meshes.emplace_back(m_objects[i]->getMesh());
    }
}
