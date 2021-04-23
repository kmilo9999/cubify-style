#include "toonmaterial.h"
#include "graphics/scene.h"
#include "graphics/camera.h"

ToonMaterial::ToonMaterial() : Material(), diffuseColor(1.0, 1.0, 1.0), specularColor(1.0, 1.0, 1.0),
    ambientColor(0, 0, 0), sssColor(0.4, 0.1, 0.1), specularPower(20.f), sssPower(.2f), sssRadius(0.1f),
    colorLevels(-0.3f, 0.4f, 1.f), colorPowers(0.1f, 0.5f, 0.9f, 1.f), colorRadius(0.05f, 0.05f, 0.02f)
{
    getShaders(this);

    // set texture
    updateProperty();
}

ToonMaterial::~ToonMaterial()
{}

void ToonMaterial::loadDiffuseTexture(const std::string &texPath)
{
    diffuseTex.load(texPath);
}

void ToonMaterial::updateProperty()
{
    m_shader->bind();
    updatePropertyWhenBind();
    m_shader->unbind();
}

void ToonMaterial::updatePropertyWhenBind()
{
    m_shader->setUniform("diffuseTex", (int)0);
    m_shader->setUniform("diffuseColor", diffuseColor);
    m_shader->setUniform("specularColor", specularColor);
    m_shader->setUniform("ambiendColor", ambientColor);
    m_shader->setUniform("sssColor", sssColor);
    m_shader->setUniform("specularPower", specularPower);
    m_shader->setUniform("sssPower", sssPower);
    m_shader->setUniform("sssRadius", sssRadius);
    m_shader->setUniform("colorLevels", colorLevels);
    m_shader->setUniform("colorPowers", colorPowers);
    m_shader->setUniform("colorRadius", colorRadius);
}

void ToonMaterial::getShaders(ToonMaterial *object)
{
    static std::shared_ptr<Shader> shader1 = std::make_shared<Shader>(":/shaders/toon.vert", ":/shaders/toon.frag");
    static std::shared_ptr<Shader> shader2 = std::make_shared<Shader>(":/shaders/outline.vert", ":/shaders/outline.frag");
    object->m_shader = shader1;
    object->m_outlineShader = shader2;
}

int ToonMaterial::getNumberOfPass() const
{
    return 2;
}

void ToonMaterial::bind(int pass)
{
    switch (pass)
    {
    case 0:
        glCullFace(GL_BACK);
        m_shader->bind();
        updatePropertyWhenBind();
        m_shader->setUniform("vp", Scene::mainScene->getCam()->getVP());
        diffuseTex.activate(0);
        break;
    case 1:
        glCullFace(GL_FRONT);
        m_outlineShader->bind();
        m_outlineShader->setUniform("vp", Scene::mainScene->getCam()->getVP());
        break;
    }
}

void ToonMaterial::unbind(int pass)
{
    switch (pass)
    {
    case 0:
        m_shader->unbind();
        break;
    case 1:
        glCullFace(GL_BACK);
        m_outlineShader->unbind();
        break;
    }
}

void ToonMaterial::setModelMatrix(const Eigen::Matrix4f &model)
{
    m_shader->bind();
    m_shader->setUniform("m", model);
    m_outlineShader->bind();
    m_outlineShader->setUniform("m", model);
    m_outlineShader->unbind();
}
