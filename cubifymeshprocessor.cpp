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
#include <igl/rotation_matrix_from_directions.h>

CubifyMeshProcessor::CubifyMeshProcessor()
{

}

void CubifyMeshProcessor::init(std::string filename)
{

    Eigen::MatrixXd V;

    Eigen::MatrixXi F;

    igl::read_triangle_mesh("./meshes/Cube3.obj", V, F);


    _shape = std::make_shared<Shape>();

      checkError();

//      Eigen::Matrix3d RM;
//      Eigen::Vector3d axis1(1,1,0);
//      Eigen::Vector3d axis2 = Eigen::Vector3d::UnitY();
//      RM = igl::rotation_matrix_from_directions(axis1, axis2);

      Eigen::MatrixXd U = V ;
//      Eigen::VectorXd _v = V.row(5);
//      Eigen::Matrix3d m;
//      m = Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitX())
//          * Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitY());
//      V.row(5) = m * _v;
      std::vector<Eigen::Matrix3d> RotationsXVertex(V.rows());

      localStep(V,U,F,RotationsXVertex);


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

void CubifyMeshProcessor::localStep(const Eigen::MatrixXd& vertices,const Eigen::MatrixXd& deform,const Eigen::MatrixXi& faces,
                                     std::vector<Eigen::Matrix3d>& rotationXvertex)
{



    Eigen::SparseMatrix<double> cotangentW;
    igl::cotmatrix(vertices,faces,cotangentW);
    Eigen::MatrixXd normals;
    igl::per_vertex_normals(vertices,faces,normals);


    std::vector<std::vector<int>> incidentFaces;
    std::vector<std::vector<int>> vertexIndicesIncidents;
    igl::vertex_triangle_adjacency(vertices.rows(),faces,incidentFaces,vertexIndicesIncidents);

    std::vector<Eigen::MatrixXi> vertexList(vertices.rows());

    std::vector<Eigen::VectorXd> weigthsVecList(vertices.rows());
    std::vector<Eigen::MatrixXd> dv(vertices.rows());





     //the barycentric area around vertex
     Eigen::SparseMatrix<double> tmp;
     igl::massmatrix(vertices,faces,igl::MASSMATRIX_TYPE_BARYCENTRIC,tmp);
     Eigen::MatrixXd areaXVertex = tmp.diagonal();

     Eigen::VectorXi b;

     std::vector<int> adjacentFaces;
     for(int i = 0 ; i <  vertices.rows(); i++)
     {
      adjacentFaces = incidentFaces[i];
      int ns = adjacentFaces.size() * 3;
      vertexList[i].resize(ns, 2);
      weigthsVecList[i].resize(ns);

      for (size_t j=0; j<adjacentFaces.size(); j++)
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

      Eigen::MatrixXi m =vertexList[i].col(0);
//      std::cout << "Here is the matrix m:" << std::endl   << m << std::endl;


      igl::slice(vertices,m,1,hV0);
      igl::slice(vertices,vertexList[i].col(1),1,hV1);
      dv[i] = (hV1 - hV0).transpose();
     }

     // Local step
    Eigen::MatrixXd z(3,vertices.rows());
    z.setZero();
    Eigen::MatrixXd u (3,vertices.rows());
    u.setZero();
    Eigen::VectorXd rhos(vertices.rows()) ;
    rhos.setConstant(1e-4);


    double lambda = 0.025;


    for(int j=0; j < vertices.rows();j++)
    //for(int j=20; j < 21;j++)
    {
        Eigen::VectorXd vz = z.col(j);
        Eigen::VectorXd uk = u.col(j);
        Eigen::VectorXd vn = normals.row(j).transpose();

        double pk = rhos(j);

        Eigen::MatrixXi hE = vertexList[j];

        Eigen::MatrixXd dU(3,hE.rows());

        Eigen::MatrixXd U_hE0, U_hE1;
        igl::slice(deform,hE.col(0),1,U_hE0);
        igl::slice(deform,hE.col(1),1,U_hE1);

        dU = (U_hE1 - U_hE0).transpose();

        Eigen::MatrixXd dV = dv[j];
        Eigen::VectorXd vertexWeigth = weigthsVecList[j];


        for (int k=0; k<100; k++)
        {

            Eigen::MatrixXd displacement = (vz-uk).transpose();
            Eigen::Matrix3d Mi;
            optimalRotationMatrix(dV,vn,pk,vertexWeigth,dU,displacement,Mi);
            Eigen::VectorXd zOld = vz;

            // New Rotation Matrix = svd3x3
            Eigen::Matrix3d matrixU;
            Eigen::Matrix<double, 3, 1> singularValues;
            Eigen::Matrix3d matrixV;
            igl::svd3x3(Mi, matrixU, singularValues, matrixV);

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
            double kIndex = lambda *  areaXVertex(j) / pk;
            Eigen::VectorXd a = newRotationM*vn+uk;

            Eigen::VectorXd lh = a.array() - kIndex;
            Eigen::VectorXd max1 = lh.array().max(0.0);

            Eigen::VectorXd rh = -a.array() - kIndex;
            Eigen::VectorXd max2 = rh.array().max(0.0);


            Eigen::VectorXd zNew = max1 - max2;

            //sum of the residuals
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf  page 15
            // ˜uk+1 ← u^k + Ri^k+1 * ~ni − z^k+1
            uk = uk + newRotationM*vn-zNew;


            // primal and dual residuals  page 34
            // subject to z − Rinˆi = 0.
            //As with all problems where the constraint is x − z = 0, the primal
            //and dual residuals take the simple form
            //r^k = x^k − z^k, s^k = −ρ(z^k − z^k−1).

            double primalResidual = (zNew - newRotationM*vn).norm();
            double dualResidual = (-pk*(zNew - zOld)).norm();


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
                pk = pk * tIncr;
                //the scaled dual variable uk = (1/ρ)yk must also be rescaled
                //after updating ρ; for example, if ρ is halved, uk should be doubled
                //before proceeding.
                uk = uk / tIncr;
            }
            else if(dualResidual > mu * primalResidual)
            {
                pk = pk / tDecr;
                //the scaled dual variable uk = (1/ρ)yk must also be rescaled
                //after updating ρ; for example, if ρ is halved, uk should be doubled
                //before proceeding.
                uk = uk * tDecr;
            }else
            {
                pk = pk;
            }

            std::cout << "RotationM" <<std::endl << newRotationM << std::endl;

            Eigen::Vector3d xx = vertices.row(j) ;
            Eigen::Vector3d newPos =  newRotationM * xx ;

           std::cout << "new pos" <<std::endl << newPos << std::endl;

            Eigen::MatrixXd m = newRotationM*dV-dU;


            // stop the optimization
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
            // page 19

            double erel = 1e-3;
            double eabs = 1e-5;

            //end condition
            double ePrim =  sqrt(2.0*(double)z.size()) * eabs + erel * std::max((newRotationM*vn).norm(),zNew.norm());
            double eDual =  sqrt((double)z.size()) * eabs + erel * (pk*uk).norm();
            if(primalResidual < ePrim &&  dualResidual < eDual)
            {
                z.col(j) = zNew;
                u.col(j) = uk;
                rhos(j) = pk;
                rotationXvertex[j] = newRotationM;
                break;

            }


        }

    }



}

void CubifyMeshProcessor::optimalRotationMatrix(const Eigen::MatrixXd &dvi,
                                                const Eigen::VectorXd &normali,
                                                double pk,
                                                const Eigen::MatrixXd &weigth,
                                                const Eigen::MatrixXd &du, const Eigen::MatrixXd &displacement,
                                                Eigen::Matrix3d& out)
{
  Eigen::Matrix3d deforfmedVertex = dvi * weigth.asDiagonal() * du.transpose();
  out = deforfmedVertex + (pk * normali * displacement);
}
