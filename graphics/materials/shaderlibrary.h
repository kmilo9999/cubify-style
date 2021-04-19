#ifndef SHADERLIBRARY_H
#define SHADERLIBRARY_H

/**
 * @brief Class to manage all shaders in the program.
 */

#include <string>
#include <map>
#include <set>
#include <GL/glew.h>

enum ShaderType
{
    VERTEX_SHADER,
    FRAGMENT_SHADER,
    GEOMETRY_SHADER,
    COMPUTE_SHADER
};

struct CachedShader
{
    std::string name;           // shader name
    ShaderType type;
    GLuint shaderid;

    unsigned int reference = 0;
};

struct CachedProgram
{
    GLuint programid;
    GLuint vertex_shader;
    GLuint geometry_shader;
    GLuint fragment_shader;
    GLuint compute_shader;

    unsigned int reference;

    CachedProgram();

    bool operator<(const CachedProgram& p2) const
    {
        if (vertex_shader != p2.vertex_shader)
            return vertex_shader < p2.vertex_shader;
        else if (fragment_shader != p2.fragment_shader)
            return fragment_shader < p2.fragment_shader;
        else if (geometry_shader != p2.geometry_shader)
            return geometry_shader < p2.geometry_shader;
        return compute_shader < p2.compute_shader;
    }
};

class ShaderLibrary
{
public:
    ~ShaderLibrary();
    ShaderLibrary(const ShaderLibrary&) = delete;
    ShaderLibrary& operator=(const ShaderLibrary&) = delete;
    static ShaderLibrary* getInstance();    // singleton

    GLuint tryGetShader(const std::string& file_name, ShaderType type);
    GLuint tryGetProgram(const std::string& vs_path, const std::string& ps_path);

    void releaseProgram(GLuint program_id);
private:
    ShaderLibrary();

    void createShaderFromSource(const std::string &source, CachedShader* item);
    void createProgramById(GLuint vs, GLuint ps, CachedProgram *item);
    std::string getFileContents(const std::string &path);

    std::map<std::string, CachedShader> m_shaders;
    std::set<CachedProgram> m_programs;
};

#endif // SHADERLIBRARY_H
