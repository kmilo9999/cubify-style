#include <Eigen/Core>       // for cross operator
#include <limits>
#include <iostream>
#include "light.h"

Light::Light()
{

}

/**
 * Directional Light
 *
 */
DirectionalLight::DirectionalLight() : m_color(Eigen::Vector3f::Zero()),
    m_direction(0, -1, 0)
{

}

DirectionalLight::DirectionalLight(const Eigen::Vector3f& color, const Eigen::Vector3f& direction) : Light(),
    m_color(color), m_direction(direction)
{
    m_direction.normalize();
}

Eigen::Vector3f DirectionalLight::getLightColor() const
{
    return m_color;
}

Eigen::Vector3f DirectionalLight::getLightDirection() const
{
    return m_direction;
}

void DirectionalLight::getLightSpaceTransform(const FrustumPoints &clippedFrustum, Eigen::Matrix4f &lightVP) const
{
    // compute light transform matrix
    // map frustum points onto the light view plane
    Eigen::Vector3f lightDir = getLightDirection();
    Eigen::Vector3f up(0, 1.f, 0);
    // in case the light is pointing along y axis
    if (std::abs(lightDir.dot(up)) > 0.95) {
        up = Eigen::Vector3f(0, 0, 1.f);
    }
    Eigen::Vector3f right = lightDir.cross(up);
    up = right.cross(lightDir);

    Eigen::Vector4f boundingBox; // min x, min y, max x, max y
    float minFloat = std::numeric_limits<float>::min();
    float maxFloat = std::numeric_limits<float>::max();
    boundingBox << maxFloat, maxFloat, minFloat, minFloat;

    float startPoint = maxFloat;
    float mean_u = 0, mean_v = 0;
    int n_rows = clippedFrustum.rows();
    for (int i = 0; i < n_rows; ++i) {
        float u = clippedFrustum.row(i).dot(right);
        float v = clippedFrustum.row(i).dot(up);

        mean_u += u;
        mean_v += v;

        if (u < boundingBox[0])
            boundingBox[0] = u;
        if (u > boundingBox[2])
            boundingBox[2] = u;
        if (v < boundingBox[1])
            boundingBox[1] = v;
        if (v > boundingBox[3])
            boundingBox[3] = v;

        float w = clippedFrustum.row(i).dot(lightDir);
        if (w < startPoint)
            startPoint = w;
    }

    mean_u /= n_rows;
    mean_v /= n_rows;

    // with this bounding range, create light view and projection matrix
    Eigen::Matrix4f lightView;
    lightView.col(0) << right, 0.f;
    lightView.col(1) << up, 0.f;
    lightView.col(2) << -lightDir, 0.f;         // because pointing to negative z
    lightView.col(3) << -mean_u, -mean_v, startPoint + 0.2f, 1.f;

    // try an example
    Eigen::Vector4f example;
    example << clippedFrustum.row(0).transpose(), 1.f;
    std::cout << lightView * example << std::endl;
}
