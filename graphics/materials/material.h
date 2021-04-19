/**
  Class for materials
  **/

#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <memory>
#include "graphics/Shader.h"

constexpr char UNIFORM_MODEL_MATRIX_NAME[] = "m";
constexpr char UNIFORM_CAMERA_MATRIX_NAME[] = "vp";

class Material
{
public:
    Material();
    virtual ~Material();

    virtual int getNumberOfPass() const {return 1;};
    virtual void bind(int pass = 0) = 0;
    virtual void unbind(int pass = 0) = 0;
    virtual void setModelMatrix(const Eigen::Matrix4f& model) = 0;
};

#endif // MATERIAL_H
