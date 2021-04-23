#ifndef TOONMATERIAL_H
#define TOONMATERIAL_H

/**
 * @brief Toon Material
 */
#include "graphics/textures/texture2d.h"
#include "material.h"
class ToonMaterial : public Material
{
public:
    ToonMaterial();
    virtual ~ToonMaterial();

    int getNumberOfPass() const override;
    void bind(int pass = 0) override;
    void unbind(int pass = 0) override;
    virtual void setModelMatrix(const Eigen::Matrix4f &model) override;

    void loadDiffuseTexture(const std::string& texPath);
    void updateProperty();

    // material property
    Eigen::Vector3f diffuseColor;
    Eigen::Vector3f specularColor;
    Eigen::Vector3f ambientColor;
    Eigen::Vector3f sssColor;

    float specularPower;
    float sssPower;
    float sssRadius;

    // some special property
    Eigen::Vector3f colorLevels;
    Eigen::Vector4f colorPowers;
    Eigen::Vector3f colorRadius;

protected:
    void createShader();
    std::shared_ptr<Shader> m_shader;
    std::shared_ptr<Shader> m_outlineShader;

    void updatePropertyWhenBind();

    Texture2D diffuseTex;

private:
    static void getShaders(ToonMaterial* object);
};

#endif // TOONMATERIAL_H
