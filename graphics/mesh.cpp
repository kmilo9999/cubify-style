#include <algorithm>
#include <igl/pathinfo.h>
#include <igl/readMESH.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>
#include "mesh.h"

std::shared_ptr<Mesh> Mesh::createQuad()
{
    // create a quad mesh
    std::shared_ptr<Mesh> ret = std::make_shared<Mesh>();
    MatrixNr v(4, 3);
    MatrixNi f(2, 3);
    v <<    -1, -1, -1,     // 3 --- 2
            1, -1, -1,      // |     |
            1, 1, -1,       // |     |
            -1, 1, -1;      // 0 --- 1
    f << 0, 1, 2,
         2, 3, 0;
    ret->update(v, f);
    return ret;
}

Mesh::Mesh() : m_vbo(0), m_ibo(0), m_vao(0), m_debugNbo(0), m_isDirty(false)
{

}

Mesh::~Mesh()
{
    GL_DELETE_BUFFER(m_vbo);
    GL_DELETE_BUFFER(m_ibo);
    GL_DELETE_BUFFER(m_debugNbo);
    GL_DELETE_VAO(m_vao);
}

void Mesh::update(MatrixNr &v, MatrixNi &f)
{
    m_v = v;
    m_f = f;
    igl::per_vertex_normals(m_v, m_f, m_norm);
    updateOpenGL();
    updateBBox();
    m_isDirty = true;
}

bool Mesh::getDirtyFlag()
{
    if (m_isDirty) {
        m_isDirty = false;
        return true;
    }
    return m_isDirty;
}

void Mesh::updateBBox()
{
    // if it's a valid mesh
    if (m_v.size() > 0 && m_f.size() > 0) {
        Eigen::Vector3f minp, maxp;
        for (size_t i = 0; i < 3; ++i) {
            minp[i] = m_v.col(i).minCoeff();
            maxp[i] = m_v.col(i).maxCoeff();
        }
        m_bbox.setBBox(minp, maxp);
    }
}

AABox Mesh::getBBox() const
{
    return m_bbox;
}

void Mesh::updateOpenGL()
{
    if (m_v.size() <= 0 || m_f.size() <= 0) {
        return;
    }
    assert (m_v.rows() == m_norm.rows());

    bool firstCreate = false;
    if (m_vbo <= 0) {
        glGenBuffers(1, &m_vbo);
        firstCreate = true;
    }

    // data type check
    MatrixNf f_v = m_v.cast<float>();
    // MatrixNu F_vbo = m_f.cast<unsigned int>();
    MatrixNf f_norm = m_norm.cast<float>();

    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    // this implementation might be slow since it asks for new memory.
    // But it's supposed to be efficient because the driver implementation should optimize it
    // And this way turns out to be switching two VBOs.
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (f_v.size() + f_norm.size() + m_uv.size()), NULL, GL_DYNAMIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * f_v.size(), f_v.data());
    glBufferSubData(GL_ARRAY_BUFFER, sizeof(float) * f_v.size(), sizeof(float) * f_norm.size(), f_norm.data());

    // since cubify doesn't change UV, so only set UV in the first call
    if (firstCreate && m_uv.size() > 0) {
        glBufferSubData(GL_ARRAY_BUFFER, sizeof(float) * (f_v.size() + f_norm.size()), sizeof(float) * m_uv.size(), m_uv.data());
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // hard coded now, in our case, ibo won't change
    if (m_ibo <= 0) {
        glGenBuffers(1, &m_ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * m_f.size(), static_cast<GLvoid*>(m_f.data()), GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

    // set up vao
    if (m_vao <= 0) {
        // Create and bind only once.
        glGenVertexArrays(1, &m_vao);
        glBindVertexArray(m_vao);
        glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, f_v.cols(), GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, f_norm.cols(), GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid*>(sizeof(float) * f_v.size()));
        if (m_uv.size() > 0) {
            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, m_uv.cols(), GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid*>(sizeof(float) * (f_v.size() + f_norm.size())));
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
}

void Mesh::draw() const
{
    if (m_vao <= 0) {
        return;
    }
    glBindVertexArray(m_vao);
    glDrawElements(GL_TRIANGLES, m_f.size(), GL_UNSIGNED_INT, reinterpret_cast<GLvoid*>(NULL));
    glBindVertexArray(0);
}

bool Mesh::load_mesh(const std::string &file_path)
{
    std::string dir, basename, ext, filename;
    igl::pathinfo(file_path, dir, basename, ext, filename);
    // Convert extension to lower case
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return std::tolower(c); });

    // now start to read
    FILE *fp = fopen(file_path.c_str(), "rb");
    if (NULL == fp)
    {
        std::cerr << "IOError:Cannot open file:" << file_path << std::endl;
        return false;
    }
    using namespace std;
    vector<vector<double > > vV,vN,vTC,vC;
    vector<vector<int > > vF,vFTC,vFN;
    vector<tuple<string, int, int>> FM;

    if (ext == "obj")
    {
        if (!igl::readOBJ(fp, vV, vTC, vN, vF, vFTC, vFN, FM))
        {
            std::cerr << "IOError:Failed to load OBJ:" << file_path << std::endl;
            fclose(fp);
            return false;
        }
        // Annoyingly obj can store 4 coordinates, truncate to xyz for this generic
        // read_triangle_mesh
        for(auto & v : vV)
        {
          v.resize(std::min(v.size(),(size_t)3));
        }
    }
    else
    {
        std::cerr << "IOError:Only obj file is supported: " << file_path << std::endl;
        fclose(fp);
        return false;
    }
    fclose(fp);

    if (vV.size() > 0)
    {
        if (!igl::list_to_matrix(vV, m_v))
        {
            std::cerr << "Cannot convert vertex array to matrix" << std::endl;
            return false;
        }
    }

    Eigen::VectorXi I,C;
    igl::polygon_corners(vF,I,C);
    Eigen::VectorXi J;
    igl::polygons_to_triangles(I,C,m_f,J);

    igl::per_vertex_normals(m_v, m_f, m_norm);
    // get uv
    if (vTC.size() > 0) {
        m_uv = MatrixNf(m_v.rows(), 2);
        for (int ri = 0; ri < m_f.rows(); ++ri) {
            int v1 = m_f(ri, 0), v2 = m_f(ri, 1), v3 = m_f(ri, 2);
            int tc1 = vFTC[ri][0], tc2 = vFTC[ri][1], tc3 = vFTC[ri][2];
            m_uv(v1, 0) = vTC[tc1][0];
            m_uv(v1, 1) = vTC[tc1][1];
            m_uv(v2, 0) = vTC[tc2][0];
            m_uv(v2, 1) = vTC[tc2][1];
            m_uv(v3, 0) = vTC[tc3][0];
            m_uv(v3, 1) = vTC[tc3][1];
        }
    }
    updateOpenGL();
    updateBBox();
    m_isDirty = true;
    return true;
}
