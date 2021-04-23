#ifndef CUBIFYMESHPROCESSOR_H
#define CUBIFYMESHPROCESSOR_H

#include "cubify/cubifydata.h"
#include "graphics/dtype.h"
#include <string>

class CubifyMeshProcessor
{
public:
    CubifyMeshProcessor();
    void init(std::string filename);

    static void localStep(const CubifyData& data,
                   std::vector<Eigen::Matrix3d>& rotationXvertex);

    static double globalStep(const CubifyData& data, std::vector<Eigen::Matrix3d>& rots,
                    Eigen::MatrixXd& Vf);

    // static bool iteration(Eigen::MatrixXd& V, Eigen::MatrixXd& U, Eigen::MatrixXi& F);
    static bool iteration(const CubifyData& data, CubifyGraphics::MatrixNd& U);
};

#endif // CUBIFYMESHPROCESSOR_H
