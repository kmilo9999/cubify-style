#include "cubifymeshprocessor.h"
#include <vector>
#include <Eigen>
#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/arap_rhs.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/slice.h>
#include <igl/min_quad_with_fixed.h>

#include "graphics/shape.h"
#include "graphics/GraphicsDebug.h"
#include <igl/svd3x3.h>

CubifyMeshProcessor::CubifyMeshProcessor()
{

}

void CubifyMeshProcessor::init(std::string filename)
{

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    igl::read_triangle_mesh("./meshes/Cube.obj", V, F);


    _shape = std::make_shared<Shape>();

      checkError();
    _shape->init(V,F);


  checkError();


 // _shape->setModelMatrix(Eigen::Affine3f(Eigen::Scaling(0.2f, 0.2f, 0.2f)));

}

void CubifyMeshProcessor::draw(Shader *m_shader)
{
  checkError();
  _shape->draw(m_shader);
  checkError();
}

void CubifyMeshProcessor::update(float seconds)
{

}

void CubifyMeshProcessor::localStep(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces,
                                     Eigen::VectorXd& energyXvertex)
{



    Eigen::SparseMatrix<double> cotangentW;
    igl::cotmatrix(vertices,faces,cotangentW);
    Eigen::MatrixXd normals;
    igl::per_vertex_normals(vertices,faces,normals);
//    Eigen::SparseMatrix<double> Ni;
//    igl::arap_rhs(V,F,V.cols(),igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, Ni);

    std::vector<std::vector<int>> incidentFaces;
    std::vector<std::vector<int>> vertexIndicesIncidents;
    igl::vertex_triangle_adjacency(vertices.rows(),faces,incidentFaces,vertexIndicesIncidents);

    std::vector<Eigen::MatrixXi> vertexList(vertices.rows());
    std::vector<Eigen::VectorXd> weigthsVecList(vertices.rows());
    std::vector<Eigen::MatrixXd> dv(vertices.rows());
     Eigen::MatrixXd U = vertices;


     //the barycentric area around vertex
     Eigen::SparseMatrix<double> tmp;
     igl::massmatrix(vertices,faces,igl::MASSMATRIX_TYPE_BARYCENTRIC,tmp);
     Eigen::MatrixXd areaXVertex = tmp.diagonal();

     Eigen::VectorXi b;

     std::vector<int> adjacentFaces;
     for(int i = 0 ; i <  vertices.rows(); i++)
     {
      adjacentFaces = incidentFaces[i];
      vertexList[i].resize(incidentFaces.size()*3, 2);
      weigthsVecList[i].resize(incidentFaces.size()*3);

      for (size_t j=0; j<incidentFaces.size(); j++)
      {
          int v0 = faces(adjacentFaces[j],0);
          int v1 = faces(adjacentFaces[j],1);
          int v2 = faces(adjacentFaces[j],2);

          vertexList[i](3*j,0) = v0;
          vertexList[i](3*j,1) = v1;
          vertexList[i](3*j+1,0) = v1;
          vertexList[i](3*j+1,1) = v2;
          vertexList[i](3*j+2,0) = v2;
          vertexList[i](3*j+2,1) = v0;


          weigthsVecList[i](3*j) = cotangentW.coeff(v0,v1);
          weigthsVecList[i](3*j+1) = cotangentW.coeff(v1,v2);
          weigthsVecList[i](3*j+2) = cotangentW.coeff(v2,v0);
      }

      dv[i].resize(3, adjacentFaces.size()*3);
      Eigen::MatrixXd hV0, hV1;
      igl::slice(vertices,vertexList[i].col(0),1,hV0);
      igl::slice(vertices,vertexList[i].col(1),1,hV1);
      dv[i] = (hV1 - hV0).transpose();
     }

     // Local step
    Eigen::MatrixXd z(vertices.rows(),3);
    z.setZero();
    Eigen::MatrixXd u (vertices.rows(),3);
    u.setZero();
    Eigen::MatrixXd rhos(vertices.rows(),3) ;
    rhos.setConstant(1e-4);

   energyXvertex.resize(vertices.rows());
   energyXvertex.setZero();

    double lambda = 0.0;


    for(int j=0; j < vertices.rows();j++)
    {
        Eigen::VectorXd vz = z.col(j);
        Eigen::VectorXd vu = u.col(j);
        Eigen::VectorXd vn = normals.row(j).transpose();

        double p = rhos(j);

        Eigen::MatrixXi hE = vertexList[j];

        Eigen::VectorXd dU(3,hE.rows());

        Eigen::MatrixXd U_hE0, U_hE1;
        igl::slice(U,hE.col(0),1,U_hE0);
        igl::slice(U,hE.col(1),1,U_hE1);

        dU = (U_hE1 - U_hE0).transpose();

        Eigen::MatrixXd dV = dv[j];
        Eigen::VectorXd vertexWeigth = weigthsVecList[j];
        Eigen::Matrix3d deforfmedVertex = dV * vertexWeigth.asDiagonal() * dU.transpose();

        for (int k=0; k<100; k++)
        {
            Eigen::Matrix3d S = deforfmedVertex + (p * vn * (vz-vu).transpose());
            S /= S.norm();
            Eigen::VectorXd zOld = vz;

            // New Rotation Matrix = svd3x3
            Eigen::Matrix3d matrixU;
            Eigen::Matrix<double, 3, 1> singularValues;
            Eigen::Matrix3d matrixV;
            igl::svd3x3(S, matrixU, singularValues, matrixV);

            double d = (matrixV * matrixU.transpose()).determinant();
            if(d < 0)
            {
                // changing the sign of the column of Ui so that det(Ri) > 0
                matrixU.col(2) *= -1;
            }
            Eigen::Matrix3d newRotationM = matrixV  * matrixU.transpose();



            //Soft Thresholding
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf page  32
            //Sκ(a)=(a − κ)+ − (−a − κ)+
            // new Z
            double kIndex = lambda *  areaXVertex(j) / p;
            Eigen::VectorXd a = newRotationM*vn+vu;

            Eigen::VectorXd lh = a.array() - kIndex;
            Eigen::VectorXd max1 = lh.array().max(0.0);

            Eigen::VectorXd rh = -a.array() - kIndex;
            Eigen::VectorXd max2 = rh.array().max(0.0);


            Eigen::VectorXd zNew = max1 - max2;

            //new U
            // u = u + R*n - z;
            Eigen::VectorXd uNew = vu + newRotationM*vn-zNew;


            // primal and dual residuals  page 34
            // subject to z − Rinˆi = 0.
            //As with all problems where the constraint is x − z = 0, the primal
            //and dual residuals take the simple form
            //r^k = x^k − z^k, s^k = −ρ(z^k − z^k−1).

            double primalResidual = (zNew - newRotationM*vn).norm();
            double dualResidual = (-p*(zNew - zOld)).norm();


            //update p and u

            //3.4.1 Varying Penalty Parameter page  20
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
            //  ρ^k+1 :=
            // τ incrρk if ||r^k||2 > µ||s^k||2
            // ρk/τ decr if ||s^k||2 > µ||r^k||2
            // ρk otherwise,

            double tIncr = 2.0,tDecr = 2.0;
            double mu = 10.0;

            if(primalResidual > mu * dualResidual)
            {
                p = p * tIncr;
                //the scaled dual variable uk = (1/ρ)yk must also be rescaled
                //after updating ρ; for example, if ρ is halved, uk should be doubled
                //before proceeding.
                vu = vu / tIncr;
            }
            else if(dualResidual > mu * primalResidual)
            {
                p = p / tDecr;
                //the scaled dual variable uk = (1/ρ)yk must also be rescaled
                //after updating ρ; for example, if ρ is halved, uk should be doubled
                //before proceeding.
                vu = vu * tDecr;
            }else
            {
                p = p;
            }


            // stop the optimization
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
            // page 19

            double erel = 1e-3;
            double eabs = 1e-5;

            double ePrim =  sqrt(2*z.size()) * eabs + erel * std::max((newRotationM*vn).norm(),zNew.norm());
            double eDual =  sqrt(z.size()) * eabs + erel * (p*vu).norm();
            if(primalResidual < ePrim &&  dualResidual < eDual)
            {
                z.col(j) = zNew;
                u.col(j) = vu;
                rhos(j) = p;

                Eigen::VectorXd m = newRotationM*dV-dU;

                // ||X|| = Trace(X*Y*transpose(X))
                double sumDiagonal = (m * vertexWeigth.asDiagonal()*m.transpose()).trace() ;

                double arap = (0.5) * sumDiagonal;


                //https://machinelearningmastery.com/vector-norms-machine-learning/
                //The L1 norm that is calculated as the sum of the absolute values of the vector.
                //https://stackoverflow.com/questions/25340940/how-do-i-compute-the-absolute-value-of-a-vector-in-eigen
                Eigen::VectorXd vectorAbs= (newRotationM*vn).cwiseAbs();
                double l1_norm =vectorAbs.sum();

                double cubeness = lambda * areaXVertex(j) *l1_norm;

                energyXvertex(j) = arap + cubeness;
                break;

            }





        }

    }



}
