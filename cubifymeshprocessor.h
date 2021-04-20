#ifndef CUBIFYMESHPROCESSOR_H
#define CUBIFYMESHPROCESSOR_H

#include "graphics/dtype.h"
#include <string>

class CubifyMeshProcessor
{
public:
    CubifyMeshProcessor();
    void init(std::string filename);

    static void globalStep(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F,Eigen::VectorXd& energyXvertex, std::vector<Eigen::Matrix3d>& rots);

    static void localStep(const Eigen::MatrixXd& V, const Eigen::MatrixXd& U,const Eigen::MatrixXi& F,
                   std::vector<Eigen::Matrix3d>& rotationXvertex);

    static void globalStep(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F,
                    Eigen::VectorXd& energyXvertex, std::vector<Eigen::Matrix3d>& rots,
                    Eigen::MatrixXd& Vf);

    static void iteration(Eigen::MatrixXd& V, Eigen::MatrixXi& F);
};

#endif // CUBIFYMESHPROCESSOR_H
