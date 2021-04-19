#ifndef AABOX_H
#define AABOX_H

/**
 * @brief The AABox class, axis-aligned bounding box, mainly used for lisp Shadow map.
 */

#include <Eigen/Dense>

class AABox
{
public:
    AABox();
    AABox(const AABox& bbox);       // copy constructor
    AABox(const Eigen::Vector3f& min, const Eigen::Vector3f& max);

    void makeContain(const AABox& box);
    void setBBox(const Eigen::Vector3f& min, const Eigen::Vector3f& max);

    AABox operator+(const AABox& p2) const {
        AABox ret;
        ret.makeContain(*this);
        ret.makeContain(p2);
        return ret;
    }

    AABox& operator+=(const AABox& p2) {
        makeContain(p2);
        return *this;
    }

    friend AABox operator*(const Eigen::Affine3f& t, const AABox& p);

    Eigen::Vector3f getMin() const { return m_min;}
    Eigen::Vector3f getMax() const { return m_max;}
    void getVertices(std::vector<Eigen::Vector3f>& ret) const;
private:
    Eigen::Vector3f m_min;
    Eigen::Vector3f m_max;
};

AABox operator*(const Eigen::Affine3f& t, const AABox& p);

#endif // AABOX_H
