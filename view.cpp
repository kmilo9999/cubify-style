#include <iostream>
#include "view.h"

#include "viewformat.h"
#include "graphics/GraphicsDebug.h"
#include "graphics/materials/toonmaterial.h"

#include <QApplication>
#include <QKeyEvent>

#include <iostream>

#include "cubifymeshprocessor.h"
using namespace std;

View::View(QWidget *parent) : QGLWidget(ViewFormat(), parent),
    m_window(parent->parentWidget()),
    m_time(), m_timer(),
    m_forward(), m_sideways(), m_vertical(),
    m_lastX(), m_lastY(),
    m_capture(false), _pause(true)
{
    // View needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    QApplication::setOverrideCursor(Qt::ArrowCursor);

    // View needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // The game loop is implemented using a timer
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(tick()));
}

View::~View()
{
    // delete m_shader;
    if (m_thread && m_thread->isRunning()) {
        m_thread->softStop();
        m_thread->terminate();
        m_thread->wait();
    }

    if (m_thread)
        m_thread->deleteLater();
}

void View::initializeGL()
{
    checkError();
    m_renderer.init();

    m_scene = std::make_shared<DemoScene>();
    Scene::mainScene = m_scene.get();

    m_camera = m_scene->getCam();
    m_scene->getCam()->setPerspective(120, width() / static_cast<float>(height()), 0.1, 50);
    m_time.start();
    m_timer.start(1000 / 60);

    // create thread
    m_thread = new CubifyThread(this);
    connect(m_thread, &CubifyThread::updateReady, this, &View::meshUpdate);

    // set inputs
    std::vector<std::shared_ptr<Mesh>> cubify_meshes;
    m_scene->getCubifyMeshes(cubify_meshes);
    m_thread->setInput(cubify_meshes);
    m_thread->start();
}

void View::paintGL()
{
    // for each render loop, update meshes
    for (size_t i = 0; i < m_updateQueue.size(); ++i) {
        CubifyData* data = m_updateQueue[i];
        data->ptr->update(data->U, data->F);
        data->updated = false;
    }
    m_updateQueue.clear();

    m_renderer.render();
}

void View::resizeGL(int w, int h)
{
    m_renderer.Resize(w, h);
    //glViewport(0, 0, w, h);
    //m_camera.setAspect(static_cast<float>(w) / h);
}

void View::mousePressEvent(QMouseEvent *event)
{
    m_capture = true;
    m_lastX = event->x();
    m_lastY = event->y();
}

void View::mouseMoveEvent(QMouseEvent *event)
{
    int deltaX = event->x() - m_lastX;
    int deltaY = event->y() - m_lastY;

    if(m_capture) {
        if(deltaX != 0 || deltaY != 0) {
            m_camera->rotate(-deltaX * 0.01f, deltaY * 0.01f);
        }
    }
    m_lastX = event->x();
    m_lastY = event->y();
}

void View::mouseReleaseEvent(QMouseEvent *)
{
    m_capture = false;
}

void View::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->delta() * 0.1f / 120;
    m_camera->zoom(zoom);
}

void View::keyPressEvent(QKeyEvent *event)
{
    // Don't remove this -- helper code for key repeat events
    if(event->isAutoRepeat()) {
        keyRepeatEvent(event);
        return;
    }

    // Feel free to remove this
    if (event->key() == Qt::Key_Escape) QApplication::quit();

    if(event->key() == Qt::Key_C) {
        m_camera->toggleOrbit();
    }
    else if(event->key() == Qt::Key_W) {
        m_forward += 1;
    }
    else if(event->key() == Qt::Key_S) {
        m_forward -= 1;
    }
    else if(event->key() == Qt::Key_A) {
        m_sideways -= 1;
    }
    else if(event->key() == Qt::Key_D) {
        m_sideways += 1;
    }
    else if(event->key() == Qt::Key_Q) {
        m_vertical -= 1;
    }
    else if(event->key() == Qt::Key_E) {
        m_vertical += 1;
    }else if(event->key() == Qt::Key_Space) {
        _pause = !_pause;
    }
    else if (event->key() == Qt::Key_P) {
        if (m_thread->isRunning()) {
            m_thread->softStop();
        }
        if (m_thread->isFinished()) {
            m_thread->start();
        }
    }
    else if(event->key() == Qt::Key_T) {
        //m_sim.toggleWire();
    }
}

void View::keyRepeatEvent(QKeyEvent *)
{
}

void View::keyReleaseEvent(QKeyEvent *event)
{
    // Don't remove this -- helper code for key repeat events
    if(event->isAutoRepeat()) {
        return;
    }

    if(event->key() == Qt::Key_W) {
        m_forward -= 1;
    }
    else if(event->key() == Qt::Key_S) {
        m_forward += 1;
    }
    else if(event->key() == Qt::Key_A) {
        m_sideways += 1;
    }
    else if(event->key() == Qt::Key_D) {
        m_sideways -= 1;
    }
    else if(event->key() == Qt::Key_Q) {
        m_vertical += 1;
    }
    else if(event->key() == Qt::Key_E) {
        m_vertical -= 1;
    }
}

void View::tick()
{
    float seconds = m_time.restart() * 0.001f;
    m_scene->tick(seconds);
    if(!_pause)
    {
     //   m_sim.update(seconds);
//

    }

    auto look = m_camera->getLook();
    look.y() = 0;
    look.normalize();
    Eigen::Vector3f perp(-look.z(), 0, look.x());
    Eigen::Vector3f moveVec = m_forward * look + m_sideways * perp + m_vertical * Eigen::Vector3f::UnitY();
    moveVec *= seconds;
    m_camera->move(moveVec);

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}

void View::meshUpdate(CubifyData* data)
{
    std::cout << "Mesh Update is called" << std::endl;
    m_updateQueue.push_back(data);
}
