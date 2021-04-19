#include "aabox.h"

AABox::AABox() : m_min(Eigen::Vector3f::Zero()), m_max(Eigen::Vector3f::Zero())
{

}

AABox::AABox(const Eigen::Vector3f& min, const Eigen::Vector3f& max) : m_min(min), m_max(max)
{

}

AABox::AABox(const AABox& bbox) : m_min(bbox.m_min), m_max(bbox.m_max)
{

}

void AABox::makeContain(const AABox &box)
{
    m_min = m_min.cwiseMin(box.m_min);
    m_max = m_max.cwiseMax(box.m_max);
}

void AABox::setBBox(const Eigen::Vector3f &min, const Eigen::Vector3f &max)
{
    m_min = min;
    m_max = max;
}

// Get 8 vertices of the bbox
void AABox::getVertices(std::vector<Eigen::Vector3f> &ret) const
{
    const Eigen::Vector3f& m = getMin();
    const Eigen::Vector3f& M = getMax();

    ret.resize(8);
    ret[0] = m;
    ret[1] = Eigen::Vector3f(M[0], m[1], m[2]);
    ret[2] = Eigen::Vector3f(M[0], M[1], m[2]);
    ret[3] = Eigen::Vector3f(m[0], M[1], m[2]);
    ret[4] = Eigen::Vector3f(m[0], m[1], M[2]);
    ret[5] = Eigen::Vector3f(M[0], m[1], M[2]);
    ret[6] = Eigen::Vector3f(M[0], M[1], M[2]);
    ret[7] = Eigen::Vector3f(m[0], M[1], M[2]);
}

AABox operator*(const Eigen::Affine3f& t, const AABox& p)
{
    AABox ret;
    Eigen::Vector3f diag = p.m_max - p.m_min;
    Eigen::Vector3f ax(diag[0], 0, 0), ay(0, diag[1], 0), az(0, 0, diag[2]);
    Eigen::Vector3f newmin = t * p.m_min;
    ax = t.linear() * ax;
    ay = t.linear() * ay;
    az = t.linear() * az;

    ret.m_min = newmin;
    ret.m_max = newmin;
    for (int dim = 0; dim < 3; ++dim) {
        if (ax[dim] >= 0)
            ret.m_max[dim] += ax[dim];
        else
            ret.m_min[dim] += ax[dim];
        if (ay[dim] >= 0)
            ret.m_max[dim] += ay[dim];
        else
            ret.m_min[dim] += ay[dim];
        if (ax[dim] >= 0)
            ret.m_max[dim] += az[dim];
        else
            ret.m_min[dim] += az[dim];
    }
    // recalculat
    return ret;
}
