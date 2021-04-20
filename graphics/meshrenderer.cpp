#include "mesh.h"
#include "Shader.h"
#include "meshrenderer.h"
#include "mesh.h"
#include "materials/material.h"
#include <iostream>

MeshRenderer::MeshRenderer() : m_mesh(nullptr), m_mat(nullptr), m_scale(1,1,1)
{
    m_pos = Eigen::Vector3f::Zero();
    m_rot = Eigen::Quaternionf::Identity();
    m_trans = Eigen::Affine3f::Identity();
}

MeshRenderer::~MeshRenderer()
{
}

bool MeshRenderer::load_mesh(const std::string &file_path)
{
    m_mesh = std::make_shared<Mesh>();
    if (!m_mesh) {
        std::cerr << "Failed to create a point to Mesh" << std::endl;
        return false;
    }
    if (!m_mesh->load_mesh(file_path)) {
        return false;       // failed to load mesh
    }
    return true;
}

void MeshRenderer::set_mesh(const std::shared_ptr<Mesh> &mesh)
{
    m_mesh = mesh;
}

void MeshRenderer::set_material(std::shared_ptr<Material> &mat)
{
    m_mat = mat;
}

AABox MeshRenderer::get_boundingbox()
{
    if (m_mesh->getDirtyFlag()) {
        m_bbox = m_trans * m_mesh->getBBox();
    }
    return m_bbox;
}

void MeshRenderer::draw() const
{
    if (!m_mesh) {
        std::cerr << "Trying to draw a mesh without mesh data" << std::endl;
        return;
    }

    if (!m_mat) {
        std::cerr << "Trying to draw without material data" << std::endl;
        return;
    }

    const int passes = m_mat->getNumberOfPass();
    m_mat->setModelMatrix(m_trans.matrix());
    for (int i = 0; i < passes; ++i) {
        m_mat->bind(i);
        m_mesh->draw();
        m_mat->unbind(i);
    }
}

MeshRenderer& MeshRenderer::set_location(const Eigen::Vector3f &pos)
{
    m_pos = pos;
    m_trans.translation() = m_pos;
    return *this;
}

MeshRenderer& MeshRenderer::set_scale(const Eigen::Vector3f &scale)
{
    m_scale = scale;
    m_trans.scale(m_scale);
    return *this;
}

std::shared_ptr<Mesh> MeshRenderer::getMesh()
{
    return m_mesh;
}
