#include "posttonemapping.h"

PostTonemapping::PostTonemapping()
{
    m_shader = std::make_unique<Shader>(":/shaders/post.vert", ":/shaders/aces.frag");
}

void PostTonemapping::execute()
{
    m_shader->bind();
    Postprocess::draw();
    m_shader->unbind();
}
