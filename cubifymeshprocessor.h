#ifndef CUBIFYMESHPROCESSOR_H
#define CUBIFYMESHPROCESSOR_H

#include <string>
#include <memory>

class Shader;
class Shape;

class CubifyMeshProcessor
{
public:
    CubifyMeshProcessor();
    void init(std::string filename);
    void draw(Shader *m_shader);
    void update(float seconds);

private:

    std::shared_ptr<Shape> _shape;
};

#endif // CUBIFYMESHPROCESSOR_H
