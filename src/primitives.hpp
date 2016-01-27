#ifndef __PRIMITIVES_HPP__
#define __PRIMITIVES_HPP__

#include <float.h>

#include <limits>
#include <utility>
#include <vector>

#include "acceleration/aabb.hpp"
#include "geometry.hpp"
#include "vector_types.hpp"

typedef Vector2f Point2f;
typedef Vector3f Point3f;
typedef Vector2d Point2d;
typedef Vector3d Point3d;

typedef std::pair<Point2f, Point2f> Edge2f;
typedef std::pair<Point3f, Point3f> Edge3f;
typedef std::pair<Point2d, Point2d> Edge2d;
typedef std::pair<Point3d, Point3d> Edge3d;

typedef std::vector<Edge2f> Curve2f;
typedef std::vector<Edge3f> Curve3f;
typedef std::vector<Edge2d> Curve2d;
typedef std::vector<Edge3d> Curve3d;

template <class T_=float>
class Primitive {
public:
  virtual AABB<T_, 3> getAABB() const = 0;
  virtual Vector<T_, 3> nearestPoint(const Vector<T_, 3> &p) const = 0; 
};

template <class T_=float> 
class VertexP : Primitive<T_> {
public:
  VertexP() {}

  VertexP(const VertexP<T_> &other) {
    (*this) = other;
  }

  void operator=(const VertexP<T_> &other) {
    pos = other.pos;
  }

  virtual AABB<T_, 3> getAABB() const {
    AABB<T_, 3> aabb;
    aabb.grow(pos);
    return aabb;
  }

  virtual Vector<T_, 3> nearestPoint(const Vector<T_, 3> &p) const {
    return pos;
  }

  Vector<T_, 3> pos;
};

template <class T_=float>
class VertexN : public VertexP<T_> {
public:
  Vector<T_, 3> nor;
};

template <class T_=float>
class VertexT : public VertexP<T_> {
public:
  Vector<T_, 2> tex;
};

template <class T_=float>
class VertexNT : public VertexP<T_>,
		 public VertexN<T_>,
		 public VertexT<T_> {};

template <template <typename S_> class Vertex_, typename T_=float>
class Triangle_ : Primitive<T_> {
public:
  Vertex_<T_> vertices[3];

  virtual AABB<T_, 3> getAABB() const {
    AABB<T_, 3> aabb;
    aabb.grow(vertices[0].pos);
    aabb.grow(vertices[1].pos);
    aabb.grow(vertices[2].pos);
    return aabb;
  }

  /*
   * Nearest point to triangle. This is just an approximation.
   * Now it is returning the closest vertex position
   */
  virtual Vector<T_, 3> nearestPoint(const Vector<T_, 3> &p) const {
    T_ d1 = (p - vertices[0].pos).length();
    T_ d2 = (p - vertices[1].pos).length();
    T_ d3 = (p - vertices[2].pos).length();

    if(d1 > d2 && d1 > d3) return vertices[0].pos;
    if(d2 > d1 && d2 > d3) return vertices[1].pos;
    return vertices[2].pos;
  }
};

template <class T_=float>
class TriangleP : public Triangle_<VertexP, T_> {};

template <class T_=float>
class TriangleN : public Triangle_<VertexN, T_> {};

template <class T_=float>
class TriangleT : public Triangle_<VertexT, T_> {};

template <class T_=float>
class TriangleNT : public Triangle_<VertexNT, T_> {};

#endif
