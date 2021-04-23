#ifndef CUBIFYMESHPROCESSOR_H
#define CUBIFYMESHPROCESSOR_H

#include "graphics/dtype.h"
#include <string>

class CubifyMeshProcessor
{
public:
    CubifyMeshProcessor();
    void init(std::string filename);

    static void localStep(const Eigen::MatrixXd& V, const Eigen::MatrixXd& U,const Eigen::MatrixXi& F,
                   std::vector<Eigen::Matrix3d>& rotationXvertex);

    static double globalStep(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F, std::vector<Eigen::Matrix3d>& rots,
                    Eigen::MatrixXd& Vf);

    static bool iteration(const Eigen::MatrixXd& V, Eigen::MatrixXd& U, const Eigen::MatrixXi& F);
};

#endif // CUBIFYMESHPROCESSOR_H
