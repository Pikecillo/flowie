#ifndef __BUFFER_OBJECT_HPP__
#define __BUFFER_OBJECT_HPP__

#include <memory>

#include <GL/glew.h>

class BufferObject {
public:
    enum Type { Array = GL_ARRAY_BUFFER,
		ElementArray = GL_ELEMENT_ARRAY_BUFFER };

    typedef std::shared_ptr<BufferObject> Ptr;

    class Desc {
    public:
	enum Type { Float = GL_FLOAT, Int = GL_INT };

	Desc(int components, Type type, unsigned int stride, long offset=0) :
	    m_components(components), m_type(type),
	    m_stride(stride), m_offset(offset) {}

	unsigned int components() const { return m_components; }
	Type type() const { return m_type; }
	unsigned int stride() const { return m_stride; }
	long offset() const { return m_offset; }

    private:
	unsigned int m_components;
	Type m_type;
	unsigned int m_stride;
	long m_offset;
    };

public:
    BufferObject(Type type) {
	m_type = type;
	glGenBuffers(1, &m_handle);
    }

    void set_content(void *data, unsigned int size_in_bytes) {
	glBindBuffer(m_type, m_handle);
	glBufferData(m_type, size_in_bytes, data, GL_STATIC_DRAW);
    }

    GLuint handle() { return m_handle; }

private:
    GLuint m_handle;
    Type m_type;
};

#endif

