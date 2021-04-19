/**
  A class used to maintain the mesh structure for both OpenGL Rendering and cubify simulation

  **/

#ifndef MESH_H
#define MESH_H

#include <GL/glew.h>
#include <memory>
#include <string>
#include "dtype.h"
#include "datatypes/aabox.h"

using namespace CubifyGraphics;

class Mesh
{
public:
    Mesh();
    virtual ~Mesh();

    void update(MatrixNr& v, MatrixNi& f);              // update data
    bool load_mesh(const std::string& file_path);       // load mesh
    void draw() const;

    void updateBBox();
    AABox getBBox() const;

    bool getDirtyFlag();

    static std::shared_ptr<Mesh> createQuad();

private:
    MatrixNr m_v;       // vertices
    MatrixNi m_f;       // faces
    MatrixNr m_norm;    // normals
    MatrixNf m_uv;      // uv coordinates
    AABox m_bbox;       // bounding box

    void updateOpenGL();

    // Data for OpenGL usage
    GLuint m_vbo;
    GLuint m_ibo;
    GLuint m_vao;
    GLuint m_debugNbo;

    bool m_isDirty;
};

#endif // MESH_H
