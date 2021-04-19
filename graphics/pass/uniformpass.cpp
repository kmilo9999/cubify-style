#include "uniformpass.h"

UniformPass::UniformPass(const std::shared_ptr<Scene>& scene) : m_scene(scene)
{
    glGenBuffers(1, &m_ubo0);

    // set up
    glBindBuffer(GL_UNIFORM_BUFFER, m_ubo0);
    glBufferData(GL_UNIFORM_BUFFER, UNIFORM_BLOCK_TRANSFORM_SIZE, NULL, GL_DYNAMIC_DRAW);
}

UniformPass::~UniformPass()
{
    if (m_ubo0 > 0) {
        glDeleteBuffers(1, &m_ubo0);
    }
}

void UniformPass::execute()
{
}
