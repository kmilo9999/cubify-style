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
        data.U = mesh->getVerticeMatrix();
        data.ptr = mesh;
        data.updated = false;
        data.finished = false;
        data.initialized = false;
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
        int finished = 0;
        while (m_data[m_index].updated || m_data[m_index].finished) {
            if (m_data[m_index].finished) {
                finished++;
                if (finished >= m_data.size()) {
                    std::cout << "Thread finish working" << std::endl;
                    m_terminate = true;
                    break;
                }
            }
            m_index = (m_index + 1) % m_data.size();
        }

        if (m_terminate) {
            break;
        }

        if (!m_data[m_index].initialized)
            m_data[m_index].initialize();

        // get the one
        CubifyGraphics::MatrixNd U;
        m_data[m_index].finished = CubifyMeshProcessor::iteration(m_data[m_index], U);
        std::cout << "Ret:" << m_data[m_index].finished;

        // store the result
        // m_data[m_index].V = V;
        m_data[m_index].U = U;

        // updated, wait for access.
        m_data[m_index].updated = true;

        // notify the main process.
        emit updateReady(&m_data[m_index]);
    }
}
