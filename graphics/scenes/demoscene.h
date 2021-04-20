#ifndef DEMOSCENE_H
#define DEMOSCENE_H

#include <set>
#include <vector>
#include <memory>
#include "graphics/mesh.h"
#include "graphics/scene.h"

class DemoScene : public Scene
{
public:
    DemoScene();

    void getCubifyMeshes(std::vector<std::shared_ptr<Mesh>>& meshes);

private:
    std::set<std::shared_ptr<MeshRenderer>> m_noCubifySet;
};

#endif // DEMOSCENE_H
