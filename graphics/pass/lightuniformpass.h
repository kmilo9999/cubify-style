#ifndef LIGHTUNIFORMPASS_H
#define LIGHTUNIFORMPASS_H

/**
 * @brief This pass collect lighting information and set the uniform buffer
 */
#include "graphics/scene.h"
#include "renderpass.h"

// max has 4 lights since we used forward rendering.
struct LightUniformBlock
{
    Eigen::Vector4f lightColor[4];
    Eigen::Vector4f lightDirection[4];
    Eigen::Vector4f lightPos[4];
    int number;     // must be the last so that the layout is continuous
    //int padding[3]; compiler will add it.

    static constexpr int buffer_size = 4 + 16 * 4 + 16 * 4 + 16 * 4;
};

class LightUniformPass : public RenderPass
{
public:
    LightUniformPass();
    ~LightUniformPass();

    virtual void execute() override;

private:
    LightUniformBlock m_data;

    GLuint m_ubo;
};

#endif // LIGHTUNIFORMPASS_H
