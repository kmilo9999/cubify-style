#ifndef VIEW_H
#define VIEW_H

//#include "simulation.h"
#include "graphics/camera.h"
#include "graphics/Shader.h"
#include "graphics/meshrenderer.h"
#include "graphics/generalpipeline.h"
#include "graphics/scenes/demoscene.h"
#include "cubify/cubifythread.h"

#include <QGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include <memory>

#include <vector>

/**
 * This is similar to your "CS1971FrontEnd" class. Here you will receive all of the input events
 * to forward to your game.
 *
 * @brief The View class
 */

class CubifyMeshProcessor;
class View : public QGLWidget
{
    Q_OBJECT

public:
    View(QWidget *parent);
    ~View();

private:
    static const int FRAMES_TO_AVERAGE = 30;

private:
    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void keyRepeatEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);

private:
    QWidget *m_window;
    QElapsedTimer m_time;
    QTimer m_timer;

    GeneralPipeline m_renderer;
    std::shared_ptr<DemoScene> m_scene;

    std::vector<CubifyData*> m_updateQueue;

    Camera *m_camera;

    int m_forward, m_sideways, m_vertical;

    int m_lastX, m_lastY;

    bool m_capture;
    bool _pause;

    CubifyThread* m_thread;

public slots:
    void meshUpdate(CubifyData* data);

private slots:
    void tick();
};

#endif // VIEW_H

