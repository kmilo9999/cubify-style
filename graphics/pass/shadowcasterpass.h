#ifndef SHADOWCASTERPASS_H
#define SHADOWCASTERPASS_H

/**
 * @brief The Shadowcasterpass class, a pass for shadows
 * Uses Lisp Shadow Map
 */

#include "renderpass.h"
#include <memory>
#include "graphics/scene.h"

class Shadowcasterpass : public RenderPass
{
public:
    Shadowcasterpass();

    virtual void execute() override;

private:;
    bool m_debug;
};

#endif // SHADOWCASTERPASS_H
