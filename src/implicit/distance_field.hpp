/*
 * University of Houston
 * Mario Rincon-Nigro. April 2013.
 */

#ifndef __DISTANCE_FIELD_HPP__
#define __DISTANCE_FIELD_HPP__

#include <iostream>

#include "acceleration/bvh.hpp"
#include "field_types.hpp"

template<class S_, int dim_>
void make_distance_field(const Vector<S_, dim_> &center, S_ radius,
                         ScalarField<S_, dim_> &field) {
  for(int i = 0; i < field.size(); i++)
    field.at(i) = (field.coordinateAt(i) - center).length() - radius;
}

template <class Primitive_, class T_>
void distance_field(const std::vector<Primitive_> &primitives,
		    ScalarField<T_, 3> &field) {
  BVH<Primitive_, T_> bvh;

  centroid_build(bvh, primitives);

  // Make sure the field covers the BVH
  AABB<T_, 3> aabb = bvh.getRoot()->getAABB();
  assert(field.getLo() <= aabb.min &&
  	 field.getHi() >= aabb.max);

  for(int i = 0; i < field.size(); i++) {
    Vector<T_, 3> p = field.coordinateAt(i);
    field.at(i) = (bvh.nearestPoint(p) - p).length();
  }
}

#endif
