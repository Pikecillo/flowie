#include <graphics/Shader.hpp>

Shader::Shader(Type type, std::string &source) {
    m_handle = glCreateShader(type);
    
    const GLchar *src = source.c_str();
    
    glShaderSource(m_handle, 1, (const GLchar**)&src, 0);
    glCompileShader(m_handle);
}

void Shader::print_info_log() {
    int length = 0;
    int written = 0;
    char *info_log;
    
    glGetShaderiv(m_handle, GL_INFO_LOG_LENGTH, &length);
    
    if(length > 0) {
	info_log = new char[length];
	glGetShaderInfoLog(m_handle, length, &written, info_log);
	std::cout << std::string(info_log) << std::endl;
	delete [] info_log;
    }
}

ShaderProgram::ShaderProgram(const Shader::Ptr &vertex_shader,
			     const Shader::Ptr &fragment_shader) :
    m_vertex_shader(vertex_shader),
    m_fragment_shader(fragment_shader) {
    
    m_handle = glCreateProgram();
    
    glAttachShader(m_handle, m_vertex_shader->handle());
    glAttachShader(m_handle, m_fragment_shader->handle());
    
    glLinkProgram(m_handle);
}

void ShaderProgram::set_attribute(const std::string &attribute_name,
				  const BufferObject::Ptr &buffer_object,
				  const BufferObject::Desc &buffer_desc) {
    GLint location;
    location = glGetAttribLocation(m_handle, attribute_name.c_str());
    glEnableVertexAttribArray(location);
    glBindBuffer(GL_ARRAY_BUFFER, buffer_object->handle());
    glVertexAttribPointer(location, buffer_desc.components(),
			  buffer_desc.type(), GL_FALSE,
			  buffer_desc.stride(),
			  (const GLvoid *)buffer_desc.offset());
}

void ShaderProgram::set_uniform_1f(const std::string &uniform_name,
				   float v0) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform1f(location, v0);
}

void ShaderProgram::set_uniform_2f(const std::string &uniform_name,
				   float v0, float v1) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform2f(location, v0, v1);
}

void ShaderProgram::set_uniform_3f(const std::string &uniform_name,
				   float v0, float v1, float v2) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform3f(location, v0, v1, v2);
}

void ShaderProgram::set_uniform_4f(const std::string &uniform_name,
				   float v0, float v1, float v2, float v3) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform4f(location, v0, v1, v2, v3);
}

void ShaderProgram::set_uniform_1i(const std::string &uniform_name,
				   int v0) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform1i(location, v0);
}

void ShaderProgram::set_uniform_2i(const std::string &uniform_name,
				   int v0, int v1) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform2i(location, v0, v1);
}

void ShaderProgram::set_uniform_3i(const std::string &uniform_name,
				   int v0, int v1, int v2) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform3i(location, v0, v1, v2);
}

void ShaderProgram::set_uniform_4i(const std::string &uniform_name,
				   int v0, int v1, int v2, int v3) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniform4i(location, v0, v1, v2, v3);
}

void ShaderProgram::set_uniform_matrix_4fv(const std::string &uniform_name,
					   unsigned int count,
					   bool transpose,
					   const float *v0) {
    int location = glGetUniformLocation(m_handle, uniform_name.c_str());
    glUniformMatrix4fv(location, count, transpose, v0);
}

void ShaderProgram::print_info_log() const {
    int length = 0;
    int written = 0;
    char *info_log;

    glGetProgramiv(m_handle, GL_INFO_LOG_LENGTH, &length);
    
    if(length > 0) {
	info_log = new char[length];
	glGetProgramInfoLog(m_handle, length, &written, info_log);
	std::cout << std::string(info_log);
	delete [] info_log;
    }
}
