#include "lightuniformpass.h"
#include <iostream>

LightUniformPass::LightUniformPass() :
    RenderPass(), m_ubo(0)
{
    // create buffer
    glGenBuffers(1, &m_ubo);
    glBindBuffer(GL_UNIFORM_BUFFER, m_ubo);
    glBufferData(GL_UNIFORM_BUFFER, LightUniformBlock::buffer_size, NULL, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glBindBufferRange(GL_UNIFORM_BUFFER, 1, m_ubo, 0, LightUniformBlock::buffer_size);

    m_data.number = 0;
    for (size_t i = 0; i < 4; ++i) {
        m_data.lightPos[i] << 0.f, 0.f, 0.f, 1.f;
        m_data.lightColor[i].setOnes();
        m_data.lightDirection[i] << 0.f, -1.f, 0.f, 0.f;
    }
}

LightUniformPass::~LightUniformPass()
{
    if (m_ubo > 0) {
        glDeleteBuffers(1, &m_ubo);
        m_ubo = 0;
    }
}

void LightUniformPass::execute()
{
    Scene* m_scene = Scene::main;
    const std::vector<std::shared_ptr<Light>>& lights = m_scene->getLights();
    m_data.number = std::min(4, (int)lights.size());

    // get lights
    for (int i = 0; i < m_data.number; ++i) {
        m_data.lightPos[i].segment(0,3) = lights[i]->getLightPosition();
        m_data.lightColor[i].segment(0,3) = lights[i]->getLightColor();
        m_data.lightDirection[i].segment(0,3) = lights[i]->getLightDirection();
    }

    // update and bind
    glBindBuffer(GL_UNIFORM_BUFFER, m_ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, LightUniformBlock::buffer_size, &m_data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}
