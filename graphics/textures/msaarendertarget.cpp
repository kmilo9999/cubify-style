#include "msaarendertarget.h"
#include <iostream>

MSAARenderTarget::MSAARenderTarget() : m_width(0), m_height(0),
    m_fbo(0), m_depthRBO(0), m_colorRBO(0), m_normalFbo{0, 0}, m_colorTbo{0, 0}, m_activeIndex(0)
{

}

MSAARenderTarget::~MSAARenderTarget()
{
    if (m_fbo > 0) {
        glDeleteFramebuffers(1, &m_fbo);
        m_fbo = 0;
    }
    if (m_depthRBO > 0) {
        glDeleteRenderbuffers(1, &m_depthRBO);
        m_depthRBO = 0;
    }
    if (m_colorRBO > 0) {
        glDeleteRenderbuffers(1, &m_colorRBO);
        m_colorRBO = 0;
    }
    if (m_normalFbo[0] > 0) {
        glDeleteFramebuffers(2, m_normalFbo);
        m_normalFbo[0] = 0;
        m_normalFbo[1] = 0;
    }
    if (m_colorTbo[0] > 0) {
        glDeleteTextures(2, m_colorTbo);
        m_colorTbo[0] = 0;
        m_colorTbo[1] = 0;
    }
}

void MSAARenderTarget::create(int width, int height)
{
    if (width <= 0 || height <= 0) {
        std::cerr << "Trying to create negative dimension MSAARenderTarget" << std::endl;
        return;
    }

    if (m_fbo > 0) {
        std::cerr << "Try to recreate a MSAARenderTarget" << std::endl;
        return;
    }

    m_width = width;
    m_height = height;

    if (m_fbo == 0) {
        glGenFramebuffers(1, &m_fbo);
    }
    if (m_depthRBO == 0) {
        glGenRenderbuffers(1, &m_depthRBO);
    }
    if (m_colorRBO == 0) {
        glGenRenderbuffers(1, &m_colorRBO);
    }
    if (m_colorTbo[0] == 0) {
        glGenTextures(2, m_colorTbo);
        // set texture information
        glBindTexture(GL_TEXTURE_2D, m_colorTbo[0]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, m_colorTbo[1]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    if (m_normalFbo[0] == 0) {
        glGenFramebuffers(2, m_normalFbo);
        glBindFramebuffer(GL_FRAMEBUFFER, m_normalFbo[0]);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_colorTbo[0], 0);
        // set depth and stencil be none
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            std::cerr << "The normal FBO for MSAARendertarget isn't complete" << std::endl;
        }
        glBindFramebuffer(GL_FRAMEBUFFER, m_normalFbo[1]);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_colorTbo[1], 0);
        // set depth and stencil be none
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            std::cerr << "The normal FBO for MSAARendertarget isn't complete" << std::endl;
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        m_activeIndex = 0;
    }

    // set up framebuffer
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

    glBindRenderbuffer(GL_RENDERBUFFER, m_depthRBO);
    glRenderbufferStorageMultisample(GL_RENDERBUFFER, SAMPLES, GL_DEPTH24_STENCIL8, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, m_depthRBO);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, m_colorRBO);
    glRenderbufferStorageMultisample(GL_RENDERBUFFER, SAMPLES, GL_RGBA16F, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, m_colorRBO);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cerr << "Incomplete Framebuffer in MSAARenderTarget" << std::endl;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void MSAARenderTarget::resize(int width, int height)
{
    if (m_fbo == 0) {
        std::cerr << "Try to resize a non-created MSAARenderTarget" << std::endl;
        return;
    }

    if (width <= 0 || height <= 0) {
        std::cerr << "Trying to resize to a non-valid dimension of MSAARenderTarget" << std::endl;
        return;
    }

    m_width = width;
    m_height = height;

    // change render buffer size
    glBindRenderbuffer(GL_RENDERBUFFER, m_depthRBO);
    glRenderbufferStorageMultisample(GL_RENDERBUFFER, SAMPLES, GL_DEPTH24_STENCIL8, width, height);
    glBindRenderbuffer(GL_RENDERBUFFER, m_colorRBO);
    glRenderbufferStorageMultisample(GL_RENDERBUFFER, SAMPLES, GL_RGB16F, width, height);
    glBindTexture(GL_TEXTURE_2D, m_colorTbo[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    glBindTexture(GL_TEXTURE_2D, m_colorTbo[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    glBindTexture(GL_TEXTURE_2D, 0);
}

bool MSAARenderTarget::hasCreated() const
{
    return (m_fbo > 0 && m_colorRBO > 0 && m_depthRBO > 0);
}

void MSAARenderTarget::bind()
{
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);
}

void MSAARenderTarget::unbind()
{
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // blit to the normal FBO
    glBindFramebuffer(GL_READ_FRAMEBUFFER, m_fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, m_normalFbo[m_activeIndex]);
    glBlitFramebuffer(0, 0, m_width, m_height, 0, 0, m_width, m_height, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void MSAARenderTarget::startPostprocessing()
{
    // bind FBO
    int alternate_index = (m_activeIndex + 1) % 2;
    glBindFramebuffer(GL_FRAMEBUFFER, m_normalFbo[alternate_index]);
    glBindTexture(GL_TEXTURE_2D, m_colorTbo[m_activeIndex]);
    m_activeIndex = alternate_index;
}

void MSAARenderTarget::endPostprocessing()
{
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void MSAARenderTarget::blit()
{
    if (m_fbo == 0)
        return;
    glBindFramebuffer(GL_READ_FRAMEBUFFER, m_fbo);
    glBlitFramebuffer(0, 0, m_width, m_height, 0, 0, m_width, m_height, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
}

void MSAARenderTarget::drawToScreen()
{
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, m_normalFbo[m_activeIndex]);
    glBlitFramebuffer(0, 0, m_width, m_height, 0, 0, m_width, m_height, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
}
