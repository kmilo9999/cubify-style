#include "cubifydata.h"
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>

void CubifyData::initialize()
{
    igl::cotmatrix(V, F, cotMatrix);
    igl::per_vertex_normals(V, F, N);

    initialized = true;
}
