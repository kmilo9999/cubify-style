#include "cubifymeshprocessor.h"
#include <Eigen>
#include <igl/read_triangle_mesh.h>
#include "graphics/shape.h"
#include "graphics/GraphicsDebug.h"

CubifyMeshProcessor::CubifyMeshProcessor()
{

}

void CubifyMeshProcessor::init(std::string filename)
{

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    igl::read_triangle_mesh("D:/cs2240/cs2240-cubify-style/cubify-style/meshes/Cube.obj", V, F);


    _shape = std::make_shared<Shape>();

      checkError();
    _shape->init(V,F);
  checkError();


  _shape->setModelMatrix(Eigen::Affine3f(Eigen::Scaling(0.2f, 0.2f, 0.2f)));

}

void CubifyMeshProcessor::draw(Shader *m_shader)
{
  checkError();
  _shape->draw(m_shader);
  checkError();
}

void CubifyMeshProcessor::update(float seconds)
{

}
