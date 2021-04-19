#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "graphics/GraphicsDebug.h"
#include "texture2d.h"

std::unique_ptr<Texture2D> Texture2D::defaultTex(nullptr);

Texture2D::Texture2D() : m_handle(0), m_width(0), m_height(0),
    m_format(0), m_type(0), m_internalFormat(0)
{

}

Texture2D::~Texture2D()
{
    clearData();
}

void Texture2D::activateDefault(int id)
{
    if (defaultTex.get() == nullptr) {
        // create default texture, a 1x1 white texture
        defaultTex = std::make_unique<Texture2D>();
        defaultTex->create(1, 1, GL_RGB, GL_RGB, GL_UNSIGNED_BYTE);
        // bind and refresh to 0
        glBindTexture(GL_TEXTURE_2D, defaultTex->m_handle);
        unsigned char data[3] = {255, 255, 255};
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, data);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    defaultTex->activate(id);
}

void Texture2D::activate(int id)
{
    if (!isValid()) {
        activateDefault(id);
        return;
    }
    glActiveTexture(GL_TEXTURE0 + id);
    glBindTexture(GL_TEXTURE_2D, m_handle);
}

bool Texture2D::isValid() const
{
    return m_handle > 0;
}

bool Texture2D::create(int width, int height, GLenum format, GLint internal_format,
                       GLenum type)
{
    clearData();
    if (width <= 0 || height <= 0) {
        std::cerr << "Invalid width " << width << " or " << height
                  << " is used to create a texture." << std::endl;
        return false;
    }

    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "GL detects an error before creating a texture " << std::endl;
        printGLErrorCodeInEnglish(err);
    }

    glGenTextures(1, &m_handle);
    glBindTexture(GL_TEXTURE_2D, m_handle);
    glTexImage2D(GL_TEXTURE_2D, 0, internal_format, width, height, 0,
                  format, type, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glBindTexture(GL_TEXTURE_2D, 0);

    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "GL detects an error when creating a texture " << std::endl;
        printGLErrorCodeInEnglish(err);
        clearData();
        return false;
    }

    m_width = width;
    m_height = height;

    m_internalFormat = internal_format;
    m_format = format;
    m_type = type;

    return m_handle > 0;
}

bool Texture2D::load(const std::string &img_file)
{
    clearData();
    int width, height, nrChannels;
    stbi_set_flip_vertically_on_load(true);
    unsigned char *data = stbi_load(img_file.c_str(), &width, &height, &nrChannels, 0);

    glGenTextures(1, &m_handle);
    glBindTexture(GL_TEXTURE_2D, m_handle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    if (data) {
        m_type = GL_UNSIGNED_INT;
        if (nrChannels == 4) {
            m_format = GL_RGBA;
            m_internalFormat = GL_RGBA;
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        }
        else if (nrChannels == 3) {
            m_format = GL_RGB;
            m_internalFormat = GL_RGB;
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        }
        else {
            std::cerr << "Currently Texture2D only supports texture with 3 or 4 channels." << std::endl;
            stbi_image_free(data);
            glBindTexture(GL_TEXTURE_2D, 0);
            clearData();
            return false;
        }
        glGenerateMipmap(GL_TEXTURE_2D);
    }
    else {
        std::cerr << "Failed to load texture " << img_file << std::endl;
        glBindTexture(GL_TEXTURE_2D, 0);
        clearData();
        return false;
    }
    stbi_image_free(data);
    glBindTexture(GL_TEXTURE_2D, 0);
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "GL detects an error when loading texture " << img_file << std::endl;
        printGLErrorCodeInEnglish(err);
        clearData();
        return false;
    }

    m_width = width;
    m_height = height;
    return true;
}

void Texture2D::clearData()
{
    if (m_handle > 0) {
        glDeleteTextures(1, &m_handle);
    }
    m_handle = 0;
    m_width = 0;
    m_height = 0;
}
