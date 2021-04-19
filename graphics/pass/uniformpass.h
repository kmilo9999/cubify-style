#ifndef UNIFORMPASS_H
#define UNIFORMPASS_H

/**
  * Uniform Pass, used to define Uniform Blocks values in shader
  */

#define UNIFORM_BLOCK_TRANSFORM 0
// two mat4 matrices
#define UNIFORM_BLOCK_TRANSFORM_SIZE 128
#define UNIFORM_BLOCK_LIGHTING 1

#include "renderpass.h"
#include "graphics/scene.h"

class UniformPass : public RenderPass
{
public:
    UniformPass(const std::shared_ptr<Scene>& scene);
    ~UniformPass();

    virtual void execute() override;

private:
    std::shared_ptr<Scene> m_scene;

    GLuint m_ubo0;
    GLuint m_ubo1;
};

#endif // UNIFORMPASS_H
