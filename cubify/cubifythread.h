#ifndef CUBIFYTHREAD_H
#define CUBIFYTHREAD_H

/**
 * @brief A thread of cubify algorithm
 */
#include <QThread>
#include <vector>
#include "cubify/cubifydata.h"
#include "graphics/mesh.h"

class CubifyThread : public QThread
{
    Q_OBJECT
    void run() override;

public:
    CubifyThread(QObject *parent = nullptr);
    void setInput(const std::vector<std::shared_ptr<Mesh>>& meshes);
    void softStop();

private:
    std::vector<CubifyData> m_data;
    int m_index;                                // current processing mesh
    bool m_terminate;                          // if we shoud stop the thread

signals:
    void updateReady(CubifyData* data);
};

#endif // CUBIFYTHREAD_H
