#ifndef RENDERTARGET_H
#define RENDERTARGET_H

/**
 * @brief A class for FrameBuffer.
 */

#include <GL/glew.h>
#include "texture2d.h"

class RenderTarget
{
public:
    RenderTarget();
    RenderTarget(int width, int height);
    ~RenderTarget();

    void create(int width, int height, bool hasColor, bool hasDepth);

private:
    GLuint m_fbo;     // frame buffer id
    Texture2D m_color;          // offscreen color texture
    Texture2D m_depthStencil;   // depth stencil texture
    int m_width, m_height;

    void clear();
};

#endif // RENDERTARGET_H
