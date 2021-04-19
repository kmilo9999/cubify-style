/**
  Base Pass, render opaque geometry
  **/

#ifndef BASEPASS_H
#define BASEPASS_H
#include "graphics/scene.h"
#include "renderpass.h"

class BasePass : public RenderPass
{
public:
    BasePass();
    virtual ~BasePass();

    virtual void execute() override;
};

#endif // BASEPASS_H
