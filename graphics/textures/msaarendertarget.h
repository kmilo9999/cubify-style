#ifndef MSAARENDERTARGET_H
#define MSAARENDERTARGET_H

/**
  * A really special FrameBuffer for MSAA GBuffer
  */

#include <GL/glew.h>

class MSAARenderTarget
{
public:
    MSAARenderTarget();
    ~MSAARenderTarget();

    void create(int width, int height);
    void resize(int width, int height);

    void bind();
    void unbind();

    void startPostprocessing();
    void endPostprocessing();

    void drawToScreen();

    void blit();

    bool hasCreated() const;

private:
    int m_width, m_height;
    GLuint m_fbo;
    GLuint m_depthRBO;
    GLuint m_colorRBO;

    GLuint m_normalFbo[2];     // used for post processing
    GLuint m_colorTbo[2];
    int m_activeIndex;

    static constexpr int SAMPLES = 4;
};

#endif // MSAARENDERTARGET_H
