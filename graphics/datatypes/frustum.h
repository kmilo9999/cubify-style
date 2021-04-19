#ifndef FRUSTUM_H
#define FRUSTUM_H

/**
 * @brief The Frustum class, used to represent the frustum
 */

class Frustum
{
public:
    Frustum();

private:
    float m_near, m_far;        // near, far plane

};

#endif // FRUSTUM_H
