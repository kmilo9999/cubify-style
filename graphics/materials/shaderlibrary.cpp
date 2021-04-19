#include "shaderlibrary.h"

#include <iostream>
#include <QFile>
#include <QString>
#include <QTextStream>

CachedProgram::CachedProgram() : programid(0), vertex_shader(0),
    geometry_shader(0), fragment_shader(0), compute_shader(0),
    reference(0)
{}

ShaderLibrary::ShaderLibrary()
{

}

ShaderLibrary::~ShaderLibrary()
{
    // try delete all shaders
    for (auto iter : m_shaders) {
        glDeleteShader(iter.second.shaderid);
    }
    m_shaders.clear();
    // try delete all programs
    for (auto iter : m_programs) {
        glDeleteProgram(iter.programid);
    }
    m_programs.clear();
}

ShaderLibrary* ShaderLibrary::getInstance()
{
    static ShaderLibrary lib;
    return &lib;
}

GLuint ShaderLibrary::tryGetShader(const std::string &file_name, ShaderType type)
{
    // check if this shader has been cached
    auto iter = m_shaders.find(file_name);
    if (iter == m_shaders.end()) {
        // no existing cache, create and compile it
        CachedShader cs;
        cs.name = file_name;
        cs.type = type;
        std::cout << "Compiling shader:" << file_name << std::endl;
        createShaderFromSource(getFileContents(file_name), &cs);
        if (cs.shaderid > 0) {
            m_shaders.insert({file_name, cs});
        }
        return cs.shaderid;
    }
    return iter->second.shaderid;
}

GLuint ShaderLibrary::tryGetProgram(const std::string &vs_path, const std::string &ps_path)
{
    GLuint vs = tryGetShader(vs_path, ShaderType::VERTEX_SHADER);
    GLuint ps = tryGetShader(ps_path, ShaderType::FRAGMENT_SHADER);
    CachedProgram program;
    program.vertex_shader = vs;
    program.fragment_shader = ps;

    std::set<CachedProgram>::iterator iter = m_programs.find(program);
    if (iter == m_programs.end()) {
        // create program
        createProgramById(vs, ps, &program);
        if (program.programid > 0) {
            m_programs.insert(program);
        }
        return program.programid;
    }
    return (*iter).programid;
}

void ShaderLibrary::createShaderFromSource(const std::string &source, CachedShader* item)
{
    GLenum shaderType = GL_VERTEX_SHADER;
    switch (item->type) {
    case VERTEX_SHADER:
        shaderType = GL_VERTEX_SHADER;
        break;
    case FRAGMENT_SHADER:
        shaderType = GL_FRAGMENT_SHADER;
        break;
    case GEOMETRY_SHADER:
        shaderType = GL_GEOMETRY_SHADER;
        break;
    case COMPUTE_SHADER:
        shaderType = GL_COMPUTE_SHADER;
        break;
    }
    GLuint shaderHandle = glCreateShader(shaderType);
    // try to compile
    const GLchar* codeArray[] = { source.c_str()};
    glShaderSource(shaderHandle, 1, codeArray, nullptr);
    glCompileShader(shaderHandle);

    // check compilation
    GLint status;
    glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        std::cerr << "Error: Could not compile shader."  << std::endl;

        GLint maxLength = 0;
        glGetShaderiv(shaderHandle, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the null character
        std::vector<GLchar> errorLog(maxLength);
        glGetShaderInfoLog(shaderHandle, maxLength, &maxLength, &errorLog[0]);

        std::cerr << &errorLog[0] << std::endl;
        item->shaderid = 0;
    } else {
        std::cout << "Shader compiled." << std::endl;
        item->shaderid = shaderHandle;
    }
}

void ShaderLibrary::createProgramById(GLuint vs, GLuint ps, CachedProgram *item)
{
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, ps);
    glLinkProgram(prog);

    // check link
    GLint linked;
    glGetProgramiv(prog, GL_LINK_STATUS, &linked);
    if (linked == GL_FALSE) {
        std::cerr << "Shader failed to link" << std::endl;

        GLint maxLength = 0;
        glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the null character
        std::vector<GLchar> errorLog(maxLength);
        glGetProgramInfoLog(prog, maxLength, &maxLength, &errorLog[0]);

        std::cerr << &errorLog[0] << std::endl;
        item->programid = 0;
        return;
    }
    item->programid = prog;
    glDetachShader(prog, vs);
    glDetachShader(prog, ps);
}

std::string ShaderLibrary::getFileContents(const std::string &path)
{
    QString qpath = QString::fromStdString(path);
    QFile file(qpath);

    if(file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream stream(&file);
        QString contents = stream.readAll();
        file.close();
        return contents.toStdString();
    }
    return "";
}

void ShaderLibrary::releaseProgram(GLuint program_id)
{
    // not implemented yet.
}
