#ifndef CUBIFYMESHPROCESSOR_H
#define CUBIFYMESHPROCESSOR_H

#include <string>
#include <Eigen>
#include <memory>

class Shader;
class Shape;
class Mesh;

class CubifyMeshProcessor
{
public:
    CubifyMeshProcessor();
   ~CubifyMeshProcessor();
    void init(std::string filename);
    void draw(Shader *m_shader);
    void update(float seconds);
    void globalStep(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F,Eigen::VectorXd& energyXvertex, std::vector<Eigen::Matrix3d>& rots);

    void localStep(const Eigen::MatrixXd& V, const Eigen::MatrixXd& U,const Eigen::MatrixXi& F,
                   std::vector<Eigen::Matrix3d>& rotationXvertex);

    void optimalRotationMatrix(const Eigen::MatrixXd& dvi, const Eigen::VectorXd& normali,
                               double pk,
                               const Eigen::MatrixXd& weigth, const Eigen::MatrixXd& du,
                               const Eigen::MatrixXd& displacement,Eigen::Matrix3d& out);

    double globalStep(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F,
                     std::vector<Eigen::Matrix3d>& rots,
                    Eigen::MatrixXd& Vf);
    void genTestRotations(const Eigen::MatrixXd& vertices,std::vector<Eigen::Matrix3d>& rots);
//    void optimalRotationMatrix(const Eigen::MatrixXd &dvi,
//                                                    const Eigen::VectorXd &normali,
//                                                    double pk,
//                                                    const Eigen::MatrixXd &weigth,
//                                                    const Eigen::VectorXd &du, const Eigen::MatrixXd &displacement,
//                                                    Eigen::Matrix3d& out);





private:

    std::shared_ptr<Shape> _shape;
    std::unique_ptr<Mesh> _mesh;


};

#endif // CUBIFYMESHPROCESSOR_H
