#include <iostream>
#include "generalpipeline.h"

#include "graphics/GraphicsDebug.h"
#include "graphics/pass/basepass.h"
#include "graphics/pass/shadowcasterpass.h"
#include "graphics/pass/lightuniformpass.h"
#include "graphics/pass/posttonemapping.h"

GeneralPipeline::GeneralPipeline()
{
    m_GBuffer = std::make_unique<MSAARenderTarget>();
}

GeneralPipeline::~GeneralPipeline()
{

}

void GeneralPipeline::init()
{
    glewExperimental = GL_TRUE;
    if(glewInit() != GLEW_OK) {
        std::cerr << "glew initialization failed" << std::endl;
    }
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // alice blue
    glClearColor(240.0f/255.0f, 248.0f/255.0f, 255.0f/255.0f, 1);

    // add passes
    m_passes.emplace_back(std::make_unique<LightUniformPass>());
    m_passes.emplace_back(std::make_unique<Shadowcasterpass>());
    m_passes.emplace_back(std::make_unique<BasePass>());

    Postprocess::InitPostprocess();
    // m_postprocess.emplace_back(std::make_unique<PostTonemapping>());

    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "pipeline initialization error:" << std::endl;
        printGLErrorCodeInEnglish(err);
    }
}

void GeneralPipeline::Resize(int width, int height)
{
    Scene* m_scene = Scene::mainScene;
    m_width = width;
    m_height = height;
    glViewport(0, 0, width, height);
    if (m_scene && m_scene->getCam()) {
        m_scene->getCam()->setAspect(static_cast<float>(width) / height);
    }
    if (m_GBuffer->hasCreated())
        m_GBuffer->resize(width, height);
    else
        m_GBuffer->create(width, height);
}

void GeneralPipeline::render()
{
    if (Scene::mainScene == nullptr) {
        std::cout << "No scene is binded to render" << std::endl;
        return;
    }
    // We render to offscreen buffer, it should be in linear space.
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "pipeline rendering error (pre-pipeline):" << std::endl;
        printGLErrorCodeInEnglish(err);
    }

    glDisable(GL_FRAMEBUFFER_SRGB);
    m_GBuffer->bind();
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    for (size_t i = 0; i < m_passes.size(); ++i) {
        m_passes[i]->execute();
    }
    m_GBuffer->unbind();

    // Post-Processing
    if (m_postprocess.size() > 0) {
        glClear(GL_COLOR_BUFFER_BIT);
        glDisable(GL_DEPTH_TEST);
        for (size_t i = 0; i < m_postprocess.size(); ++i) {
            // switch to avoid read and write to the same buffer.
            m_GBuffer->startPostprocessing();
            m_postprocess[i]->execute();
        }
    }
    m_GBuffer->endPostprocessing();

    // Output, so apply gamma correlation.
    glEnable(GL_FRAMEBUFFER_SRGB);
    m_GBuffer->drawToScreen();

    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "pipeline rendering error (post-pipeline):" << std::endl;
        printGLErrorCodeInEnglish(err);
    }
}
