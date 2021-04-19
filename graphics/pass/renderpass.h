/**
  Base Class for render passes
  **/

#ifndef RENDERPASS_H
#define RENDERPASS_H

class RenderPass
{
public:
    RenderPass();
    virtual ~RenderPass();

    virtual void execute() = 0;

};

#endif // RENDERPASS_H
