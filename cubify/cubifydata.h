#ifndef CUBIFYDATA_H
#define CUBIFYDATA_H

/**
 * @brief Data for the cubify processor
 * I separate it from the graphics shape so that it will only be modified by cubify algorithm thread.
 * Not need to care about confliction.
 */

#include <memory>
#include <Eigen/Sparse>
#include "graphics/mesh.h"
#include "graphics/dtype.h"

struct CubifyData
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    CubifyGraphics::MatrixNd U;
    bool updated;
    bool finished;
    std::shared_ptr<Mesh> ptr;

    Eigen::SparseMatrix<double> cotMatrix;
    Eigen::MatrixXd N;

    bool initialized;
    void initialize();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // CUBIFYDATA_H
