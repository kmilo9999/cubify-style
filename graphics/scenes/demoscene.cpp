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
    addPrimitive("./meshes/sphere.obj", mat)->translate(Eigen::Vector3f(-1.5f, 0.f, 0.f)).set_scale(Eigen::Vector3f(0.5f, 0.5f, 0.5f));

    std::shared_ptr<ToonMaterial> toon_apple = std::make_shared<ToonMaterial>();
    toon_apple->diffuseColor = Eigen::Vector3f(1.0f, 0.2f, 0.2f);
    toon_apple->sssColor = Eigen::Vector3f(1.f,0,0);
    toon_apple->sssPower = 0.6f;
    mat = std::dynamic_pointer_cast<Material>(toon_apple);
    addPrimitive("./meshes/apple.obj", mat);

    std::shared_ptr<ToonMaterial> pikachu_toon = std::make_shared<ToonMaterial>();
    pikachu_toon->diffuseColor = Eigen::Vector3f(1.f, 1.f, 1.f);
    pikachu_toon->ambientColor = Eigen::Vector3f(0.3f, 0.3f, 0);
    pikachu_toon->colorPowers = Eigen::Vector4f(0.3f, 0.7f, 1.f, 1.f);
    pikachu_toon->sssColor = Eigen::Vector3f(0.3f, 0.2f, 0.f);
    pikachu_toon->sssPower = 0.3f;
    pikachu_toon->loadDiffuseTexture("./meshes/pikachu.png");
    mat = std::dynamic_pointer_cast<Material>(pikachu_toon);
    addPrimitive("./meshes/pikachu.obj", mat)->translate(Eigen::Vector3f(1.5f, 0.f, 0.f));

    // mat = std::dynamic_pointer_cast<Material>(toon_bunny);
    // addPrimitive("./meshes/bunny.obj", mat)->set_location(Eigen::Vector3f(0, 2.f, 0));

    // add lighting
    std::shared_ptr<DirectionalLight> lit = std::make_shared<DirectionalLight>(Eigen::Vector3f(0.8f, 0.8f, 0.85f), Eigen::Vector3f(-1.f, -1.f, 0.f));
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
        if (m_noCubifySet.find(m_objects[i]) == m_noCubifySet.end())
            meshes.emplace_back(m_objects[i]->getMesh());
    }
}

void DemoScene::tick(float delta_time)
{
    // if you need to modify the scene.

    // rotation demo
    for (size_t i = 0; i < m_objects.size(); ++i) {
        m_objects[i]->rotate(0, 0.15f * delta_time, 0);
    }
}
