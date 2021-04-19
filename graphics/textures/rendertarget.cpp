#include "rendertarget.h"
#include <iostream>

RenderTarget::RenderTarget() : m_width(0), m_height(0), m_fbo(0)
{

}

RenderTarget::RenderTarget(int width, int height) : m_width(width), m_height(height), m_fbo(0)
{
}

RenderTarget::~RenderTarget()
{
    clear();
}

void RenderTarget::clear()
{
    m_color.clearData();
    m_depthStencil.clearData();
    if (m_fbo > 0) {
        glDeleteFramebuffers(1, &m_fbo);
        m_fbo = 0;
    }
}

void RenderTarget::create(int width, int height, bool hasColor, bool hasDepth)
{
    if (width <= 0 || height <= 0) {
        std::cerr << "Trying to create an invalid size framebuffer" << std::endl;
        return;
    }

    clear();
    m_width = width;
    m_height = height;

    glGenFramebuffers(1, &m_fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

    if (hasColor) {
        m_color.create(width, height, GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_color.m_handle, NULL);
    }

    if (hasDepth) {
        m_depthStencil.create(width, height, GL_DEPTH_COMPONENT, GL_DEPTH_COMPONENT, GL_FLOAT);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_depthStencil.m_handle, NULL);
    }

    if (hasDepth && !hasColor) {
        // Depth only buffer
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);
    }

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cerr << "Found some errors when creating framebuffer, it's incompleted" << std::endl;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
