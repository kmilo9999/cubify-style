#include "cubifymeshprocessor.h"
#include <vector>
#include <Eigen>
#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>

#include <igl/arap_rhs.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/slice.h>
#include <igl/min_quad_with_fixed.h>

#include "graphics/shape.h"
#include "graphics/GraphicsDebug.h"
#include <igl/svd3x3.h>
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOBJ.h>
#include "Mesh.h"


CubifyMeshProcessor::CubifyMeshProcessor()
{

}

CubifyMeshProcessor::~CubifyMeshProcessor()
{

}

void CubifyMeshProcessor::init(std::string filename)
{

    Eigen::MatrixXd V;

    Eigen::MatrixXi F;


    igl::read_triangle_mesh("./meshes/Cube3.obj", V, F);
//    igl::read_triangle_mesh("./meshes/Cube.obj", V, F);
//    igl::read_triangle_mesh("./meshes/bean.obj", V, F);

    _mesh = std::make_unique<Mesh>();
    _mesh->preProcess(V,F);


    Eigen::MatrixXd U = V ;

    std::vector<Eigen::Matrix3d> RotationsXVertex(V.rows());

    localStep(V,U,F,RotationsXVertex);

    std::vector<Eigen::Matrix3d> testRots;

    genTestRotations(V,testRots);

    //globalStep(V,F,E, testRots);

    Eigen::MatrixXd Vf;
    Vf.resize(V.rows(),3);
    std::cout << "RotationsXVertex" <<std::endl;

    std::cout << RotationsXVertex[0] <<std::endl;

  globalStep(V,F, testRots, Vf);


//    std::cout << V <<std::endl;

//    Eigen::Vector3d first;
//    first << V(0) ,V(8),V(16);


//    std::cout << first <<std::endl;

//    first = first.transpose() * testRots[0];
//    std::cout << first <<std::endl;
//    V(0) = first[0];
//    V(8) = first[1];
//    V(16) = first[2];

//    std::cout << V <<std::endl;



    _shape = std::make_shared<Shape>();


     checkError();

     igl::writeOBJ("./meshes/cub3CUBY.obj",Vf,F);
    _shape->init(Vf,F);



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


void CubifyMeshProcessor::genTestRotations(const Eigen::MatrixXd& vertices,
                                     std::vector<Eigen::Matrix3d>& rots)
{
    std::cout << vertices.size() <<std::endl;

    std::cout << "gen Rots" <<std::endl;
    std::vector<Eigen::Matrix3d> rotations;

    for (size_t i = 0; i <  vertices.size()/3; i++){
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        rotations.push_back(I);
    }
    //shift one vertex 15 degrees around the y axis
    Eigen::Matrix3d fifteendegreeRotPitch;
    fifteendegreeRotPitch <<  0.96f , 0.0f , 0.25f ,
                            0.0f , 1.0f , 0.0f ,
                            -0.25f , 0.0f , 0.96f;


//    rotations[0] = fifteendegreeRotPitch;

    for (size_t i = 0; i <  rotations.size(); i++){
        std::cout << rotations[i] <<std::endl <<std::endl;

    }
    rots = rotations;




}

//void CubifyMeshProcessor::preProcess(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces)
//{

//    igl::cotmatrix(vertices,faces,_cotangentWeigths);
//    Eigen::MatrixXd normals;
//    igl::per_vertex_normals(vertices,faces,_normals);


//    igl::vertex_triangle_adjacency(vertices.rows(),faces,_incidentFaces,_vertexIndicesIncidents);


//    _vertexList.resize(vertices.rows());

//    _weigthsVecList.resize(faces.rows());
//    _dv.resize(vertices.rows());

//    //the barycentric area around vertex
//    Eigen::SparseMatrix<double> tmp;
//    igl::massmatrix(vertices,faces,igl::MASSMATRIX_TYPE_BARYCENTRIC,tmp);
//    _areaXVertex = tmp.diagonal();

//    Eigen::VectorXi b;

//    std::vector<int> adjacentFaces;
//    for(int i = 0 ; i <  vertices.rows(); i++)
//    {
//     adjacentFaces = _incidentFaces[i];
//     int ns = adjacentFaces.size() * 3;
//     _vertexList[i].resize(ns, 2);
//     _weigthsVecList[i].resize(ns);

//     for (size_t j=0; j<adjacentFaces.size(); j++)
//     {
//         int v0 = faces(adjacentFaces[j],0);
//         int v1 = faces(adjacentFaces[j],1);
//         int v2 = faces(adjacentFaces[j],2);

//         _vertexList[i](3*j,0) = v0;
//         _vertexList[i](3*j,1) = v1;
//         _vertexList[i](3*j+1,0) = v1;
//         _vertexList[i](3*j+1,1) = v2;
//         _vertexList[i](3*j+2,0) = v2;
//         _vertexList[i](3*j+2,1) = v0;


//         _weigthsVecList[i](3*j) = _cotangentWeigths.coeff(v0,v1);
//         _weigthsVecList[i](3*j+1) = _cotangentWeigths.coeff(v1,v2);
//         _weigthsVecList[i](3*j+2) = _cotangentWeigths.coeff(v2,v0);
//     }

//     _dv[i].resize(3, adjacentFaces.size()*3);
//     Eigen::MatrixXd hV0, hV1;

//     Eigen::MatrixXi m =_vertexList[i].col(0);
////      std::cout << "Here is the matrix m:" << std::endl   << m << std::endl;


//     igl::slice(vertices,m,1,hV0);
//     igl::slice(vertices,_vertexList[i].col(1),1,hV1);
//     _dv[i] = (hV1 - hV0).transpose();
//    }

//}

void CubifyMeshProcessor::globalStep(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces,
                                     std::vector<Eigen::Matrix3d>& rots,
                                     Eigen::MatrixXd& Vf)
{

    int numVertices = _mesh->_numVertices;
    int numFaces = _mesh->_numFaces;


    //get cotan matrix
//    Eigen::SparseMatrix<double> cotangentW;
//    igl::cotmatrix(vertices,faces,cotangentW);

//    cotangentW*cotangentW;


    //calc Laplace Beltrami
    std::vector<std::vector<int>> incidentFaces;
    std::vector<std::vector<int>> vertexIndicesIncidents;
    igl::vertex_triangle_adjacency(vertices.rows(),faces,incidentFaces,vertexIndicesIncidents);
    std::cout << "start global" <<std::endl;


    Eigen::MatrixXd myMat(numVertices,numVertices);

    for (size_t i = 0; i < numVertices; i++){
        for (size_t j = 0; j < numVertices; j++){

            myMat ( i,j ) = 0.0;
        }
    }

   // std::cout << M.size() <<std::endl;


    for (size_t i = 0; i < incidentFaces.size(); i++){
        std::vector<int> faceInds = incidentFaces[i];

        float vertexArea = 0;

        for (size_t j = 0; j < faceInds.size(); j++){
            int faceInd = faceInds[j];
            Eigen::Vector3i vertInds;
            vertInds << faces(faceInd) , faces(faceInd+numFaces), faces(faceInd+2*numFaces);

            Eigen::Vector3d vert0;
            vert0 << vertices(vertInds[0]), vertices(vertInds[0]+numVertices), vertices(vertInds[0]+2*numVertices);


            Eigen::Vector3d vert1;
            vert1 << vertices(vertInds[1]), vertices(vertInds[1]+numVertices), vertices(vertInds[1]+2*numVertices);

            Eigen::Vector3d vert2;
            vert2 << vertices(vertInds[2]), vertices(vertInds[2]+numVertices), vertices(vertInds[2]+2*numVertices);


            Eigen::Vector3d crossP = (vert2-vert0).cross(vert1-vert0);
            float area = crossP.norm()*0.5;

            vertexArea += area;





        }
        vertexArea /= 3;
        myMat(i,i) = vertexArea;


    }
    std::cout << "Checkpoint" <<std::endl;


    std::cout << numVertices <<std::endl;

    std::cout << myMat <<std::endl;
    //std::cout << cotangentW <<std::endl;

    Eigen::MatrixXd LTEST(numVertices,numVertices);






    //LTEST = myMat*cotangentW;
    LTEST = myMat*_mesh->_cotangentWeigths;

    std::cout << myMat.inverse() <<std::endl;
    std::cout << LTEST <<std::endl;
    std::cout << LTEST.transpose() <<std::endl;

    std::cout << LTEST+LTEST.transpose() <<std::endl;



    /* GLOBAL */
    std::vector<std::vector<int>> incidentVerts;

    igl::adjacency_list(faces, incidentVerts);

    Eigen::MatrixXd Bs(numVertices,3);


    for (size_t i = 0; i < incidentVerts.size(); i++){

        std::vector<int> curNeighborInds = incidentVerts[i];

        Eigen::Vector3d P_i;
        Eigen::Matrix3d R_i = rots[i];


        P_i =  vertices.row(i);



       Eigen::VectorXd b_i(3);
       b_i = Eigen::VectorXd::Zero(3);

        for (size_t k = 0; k < curNeighborInds.size(); k++){



            int j = curNeighborInds[k];


            Eigen::Vector3d P_j;
            Eigen::Matrix3d R_j = rots[j];

              P_j =  vertices.row(j);
     //       std::cout << P_j <<std::endl;

            float curWeight = _mesh->_cotangentWeigths.coeff(i,j);


            Eigen::Vector3d P_diff = P_i - P_j;
            Eigen::Matrix3d R_sum = R_i + R_j;



            b_i += (curWeight/2.0)*R_sum*P_diff;





        }

        Bs.row(i) =  b_i;

    }


    std::cout << "Bs-------------" <<std::endl;

    std::cout << Bs <<std::endl;
    std::cout << Bs.size() <<std::endl;



    Eigen::Vector3d  curVert;
    curVert << vertices(0), vertices(8), vertices(16);
    std::cout << curVert <<std::endl;
    Eigen::MatrixXd LTESTSUM(numVertices,numVertices);

    LTESTSUM= LTEST+LTEST.transpose();

//    LTEST += LTESTT;
    std::cout << "LTESTSUM" <<std::endl;

    std::cout << LTESTSUM <<std::endl;


    for (size_t i = 0 ; i < 3; i++){


           Eigen::MatrixXd L(numVertices, numVertices);
           L = LTESTSUM/2;



           std::cout << "B" << std::endl <<std::endl;
           std::cout << Bs.col(0) << std::endl;
           std::cout << "L" << std::endl <<std::endl;

           std::cout << L << std::endl;


           Eigen::VectorXd x = L.colPivHouseholderQr().solve(Bs.col(i));
           std::cout << "The solution is:\n" << x << std::endl;
           Vf.col(i) = -x;
    }
    std::cout << "newVerts" <<std::endl;

    std::cout << Vf <<std::endl;

    std::cout << "Original" <<std::endl;

    std::cout << vertices <<std::endl;




}


void CubifyMeshProcessor::localStep(const Eigen::MatrixXd& vertices,const Eigen::MatrixXd& Uvertices,const Eigen::MatrixXi& faces,
                                     std::vector<Eigen::Matrix3d>& rotationXvertex)
{




     // Local step
//    Eigen::MatrixXd z(3,vertices.rows());
//    z.setZero();
//    Eigen::MatrixXd u (3,vertices.rows());
//    u.setZero();

    Eigen::MatrixXd  z(3,1);
    z.setZero();
    Eigen::MatrixXd u(3,1);
    u.setZero();
    Eigen::VectorXd rhos(vertices.rows()) ;
    rhos.setConstant(1e-4);


    double lambda = 0.025;


    for(int i=0; i < vertices.rows();i++)
    //for(int j=20; j < 21;j++)
    {
        Eigen::VectorXd zk = z.col(0);
        Eigen::VectorXd uk = u.col(0);
        Eigen::VectorXd vn = _mesh->_normals.row(i).transpose();

        double pk = rhos(i);

        Eigen::MatrixXi hE = _mesh->_vertexList[i];

        Eigen::MatrixXd dU(3,hE.rows());

        Eigen::MatrixXd U_hE0, U_hE1;
        igl::slice(Uvertices,hE.col(0),1,U_hE0);
        igl::slice(Uvertices,hE.col(1),1,U_hE1);

        dU = (U_hE1 - U_hE0).transpose();

        Eigen::MatrixXd dV = _mesh->_dv[i];
        Eigen::VectorXd vertexWeigth = _mesh->_weigthsVecList[i];


        //for (int j=0; j<100; j++)
        bool optimizing = true;
        while(optimizing)
        {

            Eigen::MatrixXd displacement = (zk-uk).transpose();

            Eigen::Matrix3d deforfmedVertex = dV * vertexWeigth.asDiagonal() * dU.transpose();
            Eigen::Matrix3d Mi = deforfmedVertex + (pk * vn * displacement);
           // optimalRotationMatrix(dV,vn,pk,vertexWeigth,dU,displacement,Mi);
            Eigen::VectorXd zOld = zk;

            // New Rotation Matrix = svd3x3
            Eigen::Matrix3d matrixU;
            Eigen::Matrix<double, 3, 1> singularValues;
            Eigen::Matrix3d matrixV;
            igl::svd3x3(Mi, matrixU, singularValues, matrixV);

            double d = (matrixV * matrixU.transpose()).determinant();
            if(d < 0)
            {
                // changing the sign of the column of Ui so that det(Ri) > 0
                matrixU.col(2) = -matrixU.col(2);
            }
            Eigen::Matrix3d newRotationM = matrixV  * matrixU.transpose();



            //Soft Thresholding
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf page  32
            //Sκ(a)=(a − κ)+ − (−a − κ)+
            // new Z
            double kIndex = lambda *  _mesh->_areaXVertex(i) / pk;
            Eigen::VectorXd a = newRotationM*vn+uk;

            Eigen::VectorXd lh = a.array() - kIndex;
            Eigen::VectorXd max1 = lh.array().max(0.0);

            Eigen::VectorXd rh = -a.array() - kIndex;
            Eigen::VectorXd max2 = rh.array().max(0.0);


            zk = max1 - max2;

            //sum of the residuals
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf  page 15
            // ˜uk+1 ← u^k + Ri^k+1 * ~ni − z^k+1
            uk = uk + newRotationM*vn-zk;


            // primal and dual residuals  page 34
            // subject to z − Rinˆi = 0.
            //As with all problems where the constraint is x − z = 0, the primal
            //and dual residuals take the simple form
            //r^k = x^k − z^k, s^k = −ρ(z^k − z^k−1).

            double primalResidual = (zk - newRotationM*vn).norm();
            double dualResidual = (-pk*(zk - zOld)).norm();


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

          //  std::cout << "RotationM" <<std::endl << newRotationM << std::endl;

       //     Eigen::Vector3d xx = vertices.row(j) ;
        //    Eigen::Vector3d newPos =  newRotationM * xx ;

        //   std::cout << "new pos" <<std::endl << newPos << std::endl;

    //        Eigen::MatrixXd myM = newRotationM*dV-dU;


            // stop the optimization
            //https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
            // page 19

            double erel = 1e-3;
            double eabs = 1e-5;

            //end condition
            double ePrim =  sqrt(2.0*(double)z.size()) * eabs + erel * std::max((newRotationM*vn).norm(),zk.norm());
            double eDual =  sqrt((double)z.size()) * eabs + erel * (pk*uk).norm();
            if(primalResidual < ePrim &&  dualResidual < eDual)
            {
                optimizing = false;
//                z.col(i) = zk;
//                u.col(i) = uk;
//                rhos(i) = pk;
                rotationXvertex[i] = newRotationM;

            }


        }

    }



}

//void CubifyMeshProcessor::optimalRotationMatrix(const Eigen::MatrixXd &dvi,
//                                                const Eigen::VectorXd &normali,
//                                                double pk,
//                                                const Eigen::MatrixXd &weigth,
//                                                const Eigen::MatrixXd &du, const Eigen::MatrixXd &displacement,
//                                                Eigen::Matrix3d& out)
//{
//  Eigen::Matrix3d deforfmedVertex = dvi * weigth.asDiagonal() * du.transpose();
//  out = deforfmedVertex + (pk * normali * displacement);
//}
