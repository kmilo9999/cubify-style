#ifndef MESHRENDERER_H
#define MESHRENDERER_H

/**
  A class that combines transformation, material and mesh.
  */

#include <memory>
#include <GL/glew.h>
#include "Eigen/Dense"
#include <unsupported/Eigen/OpenGLSupport>
#include "datatypes/aabox.h"

class Mesh;
class Material;

class MeshRenderer
{
public:
    MeshRenderer();
    virtual ~MeshRenderer();

    bool load_mesh(const std::string& file_path);
    void set_mesh(const std::shared_ptr<Mesh>& mesh);
    void set_material(std::shared_ptr<Material>& mat);

    Eigen::Affine3f get_transform() const {return m_trans;}
    AABox get_boundingbox();

    void draw() const;

    MeshRenderer& set_location(const Eigen::Vector3f& pos);
    MeshRenderer& set_scale(const Eigen::Vector3f& scale);
    MeshRenderer& set_rotation(const Eigen::Quaternionf rotation);
    MeshRenderer& set_rotation(float angle_x, float angle_y, float angle_z);
    MeshRenderer& translate(const Eigen::Vector3f& translation);
    MeshRenderer& rotate(float angle_x, float angle_y, float angle_z);

    Eigen::Vector3f get_scale() const;
    Eigen::Vector3f get_location() const;
    Eigen::Quaternionf get_rotation() const;

    std::shared_ptr<Mesh> getMesh();

private:
    std::shared_ptr<Mesh> m_mesh;
    std::shared_ptr<Material> m_mat;

    Eigen::Vector3f m_pos;
    Eigen::Vector3f m_scale;
    Eigen::Quaternionf m_rot;
    Eigen::Affine3f m_trans;

    AABox m_bbox;
};

#endif // MESHRENDERER_H
