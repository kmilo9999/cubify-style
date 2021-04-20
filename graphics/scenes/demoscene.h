#ifndef DEMOSCENE_H
#define DEMOSCENE_H

#include <vector>
#include <memory>
#include "graphics/mesh.h"
#include "graphics/scene.h"

class DemoScene : public Scene
{
public:
    DemoScene();

    void getCubifyMeshes(std::vector<std::shared_ptr<Mesh>>& meshes);
};

#endif // DEMOSCENE_H
