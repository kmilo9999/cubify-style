/**
  Rendering Pipeline for the project
  A naive forward renderer with postprocessing pass.

  Support up to 4 dynamic lights and 1 main directional light.
  **/

#ifndef GENERALPIPELINE_H
#define GENERALPIPELINE_H

/**
  * This is a forward rendering pipeline which supports at maximum 4 light sources and postprocessing
  * MSAA is implemented for anti-aliasing
  */

#include <memory>
#include <vector>
#include "graphics/pass/renderpass.h"
#include "graphics/pass/postprocess.h"
#include <GL/glew.h>

#include "graphics/textures/msaarendertarget.h"     // used for off-screen buffer and further postprocessing.

class Scene;

class GeneralPipeline
{
public:
    GeneralPipeline();
    ~GeneralPipeline();

    void init();
    void Resize(int width, int height);
    void render();

protected:
    int m_width, m_height;
    std::vector<std::unique_ptr<RenderPass>> m_passes;
    std::vector<std::unique_ptr<Postprocess>> m_postprocess;
    std::unique_ptr<MSAARenderTarget> m_GBuffer;
};

#endif // GENERALPIPELINE_H
