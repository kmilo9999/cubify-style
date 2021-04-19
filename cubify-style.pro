QT += core gui opengl

TARGET = cubify-style
TEMPLATE = app

CONFIG += c++17 console
CONFIG -= app_bundle

unix:!macx {
    LIBS += -lGLU
}
win32 {
    DEFINES += GLEW_STATIC
    LIBS += -lopengl32 -lglu32
}

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        cubifymeshprocessor.cpp \
        graphics/datatypes/aabox.cpp \
        graphics/datatypes/frustum.cpp \
        graphics/generalpipeline.cpp \
        graphics/lights/light.cpp \
        graphics/materials/material.cpp \
        graphics/materials/shaderlibrary.cpp \
        graphics/materials/toonmaterial.cpp \
        graphics/mesh.cpp \
        graphics/meshrenderer.cpp \
        graphics/pass/basepass.cpp \
        graphics/pass/lightuniformpass.cpp \
        graphics/pass/postprocess.cpp \
        graphics/pass/posttonemapping.cpp \
        graphics/pass/renderpass.cpp \
        graphics/pass/shadowcasterpass.cpp \
        graphics/pass/uniformpass.cpp \
        graphics/scene.cpp \
        graphics/textures/msaarendertarget.cpp \
        graphics/textures/rendertarget.cpp \
        graphics/textures/texture2d.cpp \
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
    cubifymeshprocessor.h \
    graphics/GraphicsDebug.h \
    graphics/MeshLoader.h \
    graphics/Shader.h \
    graphics/ShaderAttribLocations.h \
    graphics/camera.h \
    graphics/datatypes/aabox.h \
    graphics/datatypes/frustum.h \
    graphics/dtype.h \
    graphics/generalpipeline.h \
    graphics/lights/light.h \
    graphics/materials/material.h \
    graphics/materials/shaderlibrary.h \
    graphics/materials/toonmaterial.h \
    graphics/mesh.h \
    graphics/meshrenderer.h \
    graphics/pass/basepass.h \
    graphics/pass/lightuniformpass.h \
    graphics/pass/postprocess.h \
    graphics/pass/posttonemapping.h \
    graphics/pass/renderpass.h \
    graphics/pass/shadowcasterpass.h \
    graphics/pass/uniformpass.h \
    graphics/scene.h \
    graphics/shape.h \
    graphics/textures/msaarendertarget.h \
    graphics/textures/rendertarget.h \
    graphics/textures/stb_image.h \
    graphics/textures/texture2d.h \
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
