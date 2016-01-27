#ifndef __SHADER_HPP__
#define __SHADER_HPP__

#include <memory>
#include <iostream>
#include <string>

#include <graphics/BufferObject.hpp>

#include <GL/glew.h>

class Shader {
public:
    enum Type { Vertex = GL_VERTEX_SHADER,
		Fragment = GL_FRAGMENT_SHADER };

    typedef std::shared_ptr<Shader> Ptr;

public:
    Shader(Type type, std::string &source);

    GLuint handle() const { return m_handle; }

    void print_info_log();

private:
    GLuint m_handle;
};

class ShaderProgram {
public:
    typedef std::shared_ptr<ShaderProgram> Ptr;

public:
    ShaderProgram(const Shader::Ptr &vertex_shader,
		  const Shader::Ptr &fragment_shader);

    GLuint handle() const { return m_handle; }

    void set_attribute(const std::string &attribute_name,
		       const BufferObject::Ptr &buffer_object,
		       const BufferObject::Desc &buffer_desc);
    void set_uniform_1f(const std::string &uniform_name,
			float v0);
    void set_uniform_2f(const std::string &uniform_name,
			float v0, float v1);
    void set_uniform_3f(const std::string &uniform_name,
			float v0, float v1, float v2);
    void set_uniform_4f(const std::string &uniform_name,
			float v0, float v1, float v2, float v3);
    void set_uniform_1i(const std::string &uniform_name,
			int v0);
    void set_uniform_2i(const std::string &uniform_name,
			int v0, int v1);
    void set_uniform_3i(const std::string &uniform_name,
			int v0, int v1, int v2);
    void set_uniform_4i(const std::string &uniform_name,
			int v0, int v1, int v2, int v3);
    void set_uniform_matrix_4fv(const std::string &uniform_name,
				unsigned int count,
				bool transpose,
				const float *v0);

    void print_info_log() const;

private:
    Shader::Ptr m_vertex_shader;
    Shader::Ptr m_fragment_shader;
    GLuint m_handle;
};

#endif
