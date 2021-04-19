#ifndef POSTTONEMAPPING_H
#define POSTTONEMAPPING_H

#include "graphics/Shader.h"
#include "postprocess.h"

/**
 * @brief Tone Mapping Postprocess
 */

class PostTonemapping : public Postprocess
{
public:
    PostTonemapping();

    virtual void execute() override;
private:
    std::unique_ptr<Shader> m_shader;
};

#endif // POSTTONEMAPPING_H
