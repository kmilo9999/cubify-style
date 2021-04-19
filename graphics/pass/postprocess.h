#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "graphics/materials/material.h"
#include "renderpass.h"

/**
 * Special render pass, which by default render a quad using the offscreen texture.
 */

class Mesh;     // used to draw quad

class Postprocess : public RenderPass
{
public:
    Postprocess();
    static void InitPostprocess();      // called by RP to set up static variables

protected:
    void draw();                        // draw the quad

private:
    static std::shared_ptr<Mesh> QuadMesh;
};

#endif // POSTPROCESS_H
