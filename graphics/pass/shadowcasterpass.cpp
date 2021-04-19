#include "shadowcasterpass.h"
#include "graphics/dtype.h"
#include "graphics/mesh.h"
#include "graphics/meshrenderer.h"
#include "graphics/materials/toonmaterial.h"
#include <iostream>
Shadowcasterpass::Shadowcasterpass() : RenderPass()
{
    m_debug = false;
}

void Shadowcasterpass::execute()
{
    /*
    if (m_debug)
        return;         // for debug, only run it once.
    m_debug = true;
    // first make sure the scene isn't empty.
    AABox scenebbox = m_scene->getBoundingBox();
    const Eigen::Vector3f min = scenebbox.getMin();
    const Eigen::Vector3f max = scenebbox.getMax();

    // check if the scene bounding box is empty
    if ((max - min).squaredNorm() <= 1e-6f) {
        // empty scene, no need to render shadows
        return;
    }

    // clip the view frustum with the scene bounding box
    // get camera
    Camera* cam = m_scene->getCam();
    if (!cam) {
        std::cerr << "Error: scene->getCam() returns nullptr" << std::endl;
        return;
    }

    // compute the eight vertices
    Eigen::Matrix4f camVP = cam->getProjection() * cam->getView();
    Eigen::Matrix4f invVP = camVP.inverse();
    Eigen::Matrix<float, 8, 4, Eigen::RowMajor> points;
    FrustumPoints frustum;
    points <<   -1, -1, -1, 1,
                1, -1, -1, 1,
                1, 1, -1, 1,
                -1, 1, -1, 1,
                -1, -1, 1, 1,
                1, -1, 1, 1,
                1, 1, 1, 1,
                -1, 1, 1, 1;

    for (size_t i = 0; i < 8; ++i) {
        Eigen::Vector4f transformed = invVP * points.row(i).transpose();
        transformed /= transformed.w();
        frustum.row(i) = transformed.head(3);
    }

    // compute clipped points
    for (size_t i = 0; i < 4; ++i) {
        // compute a ray
        Eigen::Vector3f start = frustum.row(i);
        Eigen::Vector3f end = frustum.row(i + 4);
        Eigen::Vector3f dir = end - start;
        // 3 axes clip
        float start_t = 2.f, end_t = -1.f;

        Eigen::Vector3f t1 = (min - start).array() / dir.array();
        Eigen::Vector3f t2 = (max - start).array() / dir.array();

        // get valid min and max t.
        for (int ti = 0; ti < 3; ++ti) {
            if (t1[ti] >= 0.f && t1[ti] <= 1.f) {
                start_t = std::min(start_t, t1[ti]);
                end_t = std::max(end_t, t1[ti]);
            }
            if (t2[ti] >= 0.f && t2[ti] <= 1.f) {
                start_t = std::min(start_t, t2[ti]);
                end_t = std::max(end_t, t2[ti]);
            }
        }

        if (start_t > end_t) {          // means it's inside the bbox or totally outside of it.
            start_t = 0;
            end_t = 1.f;
        }

        frustum.row(i) = start + start_t * dir;
        frustum.row(i + 4) = start + end_t * dir;
    }

    // create face index to draw the frustum
    CubifyGraphics::MatrixNi indices(12, 3);
    indices <<  0, 1, 2,
                0, 2, 3,
                1, 5, 6,
                1, 6, 2,
                4, 0, 3,
                4, 3, 7,
                4, 5, 1,
                4, 1, 0,
                3, 2, 6,
                3, 6, 7,
                5, 4, 7,
                5, 7, 6;

    std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();
    CubifyGraphics::MatrixNr verts = frustum.cast<CubifyGraphics::real>();
    mesh->update(verts, indices);

    std::shared_ptr<MeshRenderer> render = std::make_shared<MeshRenderer>();
    render->set_mesh(mesh);

    std::shared_ptr<Material> mat = std::make_shared<ToonMaterial>();
    render->set_material(mat);

    m_scene->addPrimitive(render);

    // get light sources
    const std::vector<std::shared_ptr<Light>>& lights = m_scene->getLights();
    // iterate for most 4 lights
    for (size_t li = 0; li < 4; ++li) {
        const std::shared_ptr<Light>& light = lights[li];
        if (light->isShadowCaster()) {
            Eigen::Matrix4f lightVP;
            light->getLightSpaceTransform(points, lightVP);

            // TODO: break, only support one shadow caster.
            break;
        }
    }*/
}
