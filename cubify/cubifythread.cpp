#include <iostream>
#include <set>
#include "cubifymeshprocessor.h"
#include "cubifythread.h"

CubifyThread::CubifyThread(QObject *parent) : QThread(parent)
{}

void CubifyThread::setInput(const std::vector<std::shared_ptr<Mesh> > &meshes)
{
    if (isRunning()) {
        std::cout << "Cannot set input when the thread is running." << std::endl;
        return;
    }

    m_data.clear();
    std::set<std::shared_ptr<const Mesh>> saved;
    // check if it's already in the queue
    for (size_t i = 0; i < meshes.size(); ++i)
    {
        std::shared_ptr<Mesh> mesh = meshes[i];
        if (saved.find(mesh) != saved.end())
            continue;           // already has it
        // else add it
        CubifyData data;
        data.V = mesh->getVerticeMatrix();
        data.F = mesh->getFaceMatrix();
        data.ptr = mesh;
        data.updated = false;
        m_data.push_back(data);
        saved.insert(mesh);
    }
}

void CubifyThread::softStop()
{
    if (isRunning())
        m_terminate = true;
}

void CubifyThread::run()
{
    m_index = 0;
    m_terminate = false;

    if (m_data.empty()) {
        std::cout << "CubifyThread got 0 inputs so it stops" << std::endl;
        return;
    }
    // this is the main process
    while(!m_terminate) {
        // pick out one
        while (m_data[m_index].updated) {
            m_index = (m_index + 1) % m_data.size();
        }

        // get the one
        Eigen::MatrixXd V = m_data[m_index].V;
        Eigen::MatrixXi F = m_data[m_index].F;
        CubifyMeshProcessor::iteration(V, F);

        // store the result
        m_data[m_index].V = V;

        // updated, wait for access.
        m_data[m_index].updated = true;

        // notify the main process.
        emit updateReady(&m_data[m_index]);
    }
}
