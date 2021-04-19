#include "postprocess.h"
#include "graphics/mesh.h"

std::shared_ptr<Mesh> Postprocess::QuadMesh(nullptr);

Postprocess::Postprocess()
{

}

void Postprocess::InitPostprocess()
{
    if (!QuadMesh) {
        QuadMesh = Mesh::createQuad();
    }
}

void Postprocess::draw()
{
    QuadMesh->draw();
}
