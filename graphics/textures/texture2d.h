#ifndef TEXTURE2D_H
#define TEXTURE2D_H

/**
 * @brief Class for 2D textures
 */

#include <memory>
#include <string>
#include <GL/glew.h>

class Texture2D
{
public:
    Texture2D();
    ~Texture2D();

    bool create(int width, int height, GLenum format, GLint internal_format,
                GLenum type);
    bool load(const std::string& img_file);

    void activate(int id);
    bool isValid() const;
    void clearData();

    friend class RenderTarget;

    static void activateDefault(int id);

private:
    GLuint m_handle;
    int m_width, m_height;
    GLenum m_format, m_type;
    GLint m_internalFormat;

    static std::unique_ptr<Texture2D> defaultTex;
};

#endif // TEXTURE2D_H
