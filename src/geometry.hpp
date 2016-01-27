/*
 * University of Houston
 * Mario Rincon-Nigro. April 2013.
 */

#ifndef __GEOMETRY_HPP__
#define __GEOMETRY_HPP__

#include <iostream>
#include <limits>
#include <vector>

#include "vector_types.hpp"

template <class T_>
class Line {
public:
  Line(const Vector<T_, 3> &point, const Vector<T_, 3> &vector)
    : p(point), v(vector) {}
  Vector<T_, 3> p;
  Vector<T_, 3> v;
};

template <class T_>
class Segment {
public:
  Segment(const Vector<T_, 3> &q0, const Vector<T_, 3> &q1)
    : p0(q0), p1(q1) {}

  Vector<T_, 3> p0;
  Vector<T_, 3> p1;
};

template <class T_>
class Plane {
public:
  Plane(const Vector<T_, 3> &point, const Vector<T_, 3> &normal)
    : p(point), n(normal) {}
  Vector<T_, 3> p;
  Vector<T_, 3> n;
};

/*
 * Minimum distance from point to line
 */
template <class T_>
T_ min_distance(const Line<T_> &line, const Vector<T_, 3> &p) {
  return (line.v.cross(line.p - p)).length() / line.v.length();
}

/*
 * Minimum distance from point to edge
 */
template <class T_>
T_ min_distance(const Segment<T_> &s, const Vector<T_, 3> &p) {
  T_ d = min_distance(Line<T_>(s.p0, s.p1 - s.p0), p);

  if((p - s.p0).dot(s.p1 - s.p0) < 0) d = (p - s.p0).length();
  if((p - s.p1).dot(s.p0 - s.p1) < 0) d = (p - s.p1).length();
  
  return d;
}

/*
 * Minimum distance from point to plane
 */
template <class T_>
T_ min_distance(const Plane<T_> &s, const Vector<T_, 3> &p) {
  return s.n.project(p - s.p);
}

#endif
