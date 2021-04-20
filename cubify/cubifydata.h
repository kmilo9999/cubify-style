#ifndef CUBIFYDATA_H
#define CUBIFYDATA_H

/**
 * @brief Data for the cubify processor
 * I separate it from the graphics shape so that it will only be modified by cubify algorithm thread.
 * Not need to care about confliction.
 */

#include <memory>
#include "graphics/mesh.h"
#include "graphics/dtype.h"

struct CubifyData
{
    CubifyGraphics::MatrixNd V;
    CubifyGraphics::MatrixNi F;
    bool updated;
    std::shared_ptr<Mesh> ptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // CUBIFYDATA_H
