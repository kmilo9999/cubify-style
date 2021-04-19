/**
  Header file that defines some types used for graphics
  **/

#ifndef DTYPE_H
#define DTYPE_H

#include <Eigen/Dense>

#define GL_DELETE_BUFFER(x) if(x>0) {glDeleteBuffers(1,&x); x=0;}
#define GL_DELETE_VAO(x) if(x>0) {glDeleteVertexArrays(1,&x); x=0;}

namespace CubifyGraphics
{
using float32 = float;
using float64 = double;

#ifndef FLOAT_STORE
using real = double;
static constexpr double PI = 3.14159265358979323846;
static constexpr double INV_PI = 1.0 / PI;
static constexpr double INV_2PI = INV_PI / 2.0;

static constexpr bool realIsDouble = true;

#else
using real = float;
static constexpr float PI = 3.14159265359f;
static constexpr float INV_PI = 1.f / PI;
static constexpr float INV_2PI = INV_PI / 2.f;

static constexpr bool realIsDouble = false;
#endif

// from https://github.com/hi2p-perim/lightmetrica-v2
real constexpr operator"" _f(long double v) {
    return real(v);
}

real constexpr operator"" _f(unsigned long long v) {
    return real(v);
}

using MatrixNr = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixNd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixNf = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixNi = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixNu = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
}


#endif // DTYPE_H
