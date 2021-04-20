#include "basepass.h"

BasePass::BasePass()
{

}

BasePass::~BasePass()
{

}

void BasePass::execute()
{
    // iterate all objects in the scene and render them
    Scene* m_scene = Scene::mainScene;
    auto drawables = m_scene->getPrimitives();
    for (auto primitive : drawables) {
        primitive->draw();
    }
}
