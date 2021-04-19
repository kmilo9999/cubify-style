#ifndef LIGHT_H
#define LIGHT_H

/**
 * @brief Interfaces for light
 */

#include <Eigen/Dense>
#include <unsupported/Eigen/OpenGLSupport>

#include "../datatypes/aabox.h"

using FrustumPoints = Eigen::Matrix<float, 8, 3, Eigen::RowMajor>;

class Light
{
public:
    Light();

    virtual bool isDirectional() const = 0;
    virtual bool isShadowCaster() const = 0;
    virtual Eigen::Vector3f getLightColor() const = 0;
    virtual Eigen::Vector3f getLightDirection() const = 0;
    virtual Eigen::Vector3f getLightPosition() const = 0;
    virtual void getLightSpaceTransform(const FrustumPoints& clippedFrustum, Eigen::Matrix4f& lightVP) const = 0;
};

class DirectionalLight : public Light
{
public:
    DirectionalLight();
    DirectionalLight(const Eigen::Vector3f& color, const Eigen::Vector3f& direction);

    virtual bool isDirectional() const override { return true; }
    virtual bool isShadowCaster() const override { return true; }
    virtual Eigen::Vector3f getLightColor() const override;
    virtual Eigen::Vector3f getLightPosition() const override { return Eigen::Vector3f::Zero(); }
    virtual Eigen::Vector3f getLightDirection() const override;
    virtual void getLightSpaceTransform(const FrustumPoints& clippedFrustum, Eigen::Matrix4f& lightVP) const override;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    Eigen::Vector3f m_color;
    Eigen::Vector3f m_direction;
};

#endif // LIGHT_H
