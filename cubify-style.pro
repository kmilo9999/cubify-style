QT += core gui opengl

TARGET = cubify-style
TEMPLATE = app

QMAKE_CXXFLAGS += -mstackrealign -Wa,-mbig-obj
CONFIG += c++17

unix:!macx {
    LIBS += -lGLU
}
win32 {
    DEFINES += GLEW_STATIC
    LIBS += -lopengl32 -lglu32
}
SOURCES += \
        Mesh.cpp \
        cubifymeshprocessor.cpp \
        libs/glew-1.10.0/src/glew.c \
        graphics/GraphicsDebug.cpp \
        graphics/MeshLoader.cpp \
        graphics/Shader.cpp \
        graphics/camera.cpp \
        graphics/shape.cpp \
        main.cpp \
        mainwindow.cpp \
        view.cpp \
        viewformat.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    Mesh.h \
    cubifymeshprocessor.h \
    graphics/GraphicsDebug.h \
    graphics/MeshLoader.h \
    graphics/Shader.h \
    graphics/ShaderAttribLocations.h \
    graphics/camera.h \
    graphics/shape.h \
    mainwindow.h \
    ui_mainwindow.h \
    view.h \
    viewformat.h

FORMS += \
    mainwindow.ui


RESOURCES += \
    res/shaders/shaders.qrc \
    res/shaders/shaders.qrc

DISTFILES += \
    res/shaders/shader.frag \
    res/shaders/shader.vert \
    res/shaders/shader.vert \
    res/shaders/shader.frag


INCLUDEPATH += libs glm libs/glew-1.10.0/include libs/Eigen/ libs/libigl/include libs/libigl/external
DEPENDPATH += libs glm libs/glew-1.10.0/include libs/Eigen/  libs/libigl/include libs/libigl/external
