/*
 * University of Houston
 * Mario Rincon-Nigro. May 2013.
 */

#ifndef __RECONSTRUCTION_HPP__
#define __RECONSTRUCTION_HPP__

#include <limits>

#include "field_types.hpp"
#include "acceleration/bvh.hpp"
#include "mesh/triangle_soup.hpp"
#include "numerical/rbf.hpp"

/*
 * Fits a PolygonSoupOBJ using an RBF
 */
template <class S_>
S_ rbf_surface(const PolygonSoupOBJ<S_> &obj_mesh,
                 RBF<S_, 3> &rbf) {
  std::vector<VertexP<S_> > vprimitives;
  BVH<VertexP<S_>, S_> vbvh;
  std::vector<S_> distances;
  std::vector<Vector<S_, 3> > points;

  primitive_sequence(obj_mesh, vprimitives);
  centroid_build(vbvh, vprimitives);

  std::cout << "Num centers: " << obj_mesh.positions.size() << std::endl;

  for(unsigned int i = 0; i < obj_mesh.positions.size(); i++) {
    Vector<S_, 3> p, p_in, p_out, n;

    p = obj_mesh.positions[i];
    n = obj_mesh.normals[i];

    S_ t = 1.0f;
    p_in = p - n * t;
    while(vbvh.nearestPoint(p_in) != p) {
      t *= 0.9f;
      p_in = p - n * t;
    }

    t = 1.0f;
    p_out = p + n * t;
    while(vbvh.nearestPoint(p_out) != p) {
      t *= 0.9f;
      p_out = p + n * t;
    }

    points.push_back(p);
    points.push_back(p_in);
    points.push_back(p_out);
    distances.push_back(0.0f);
    distances.push_back(-(p - p_in).length());
    distances.push_back((p - p_out).length());
  }

  // Interpolating with an RBF
  rbf = RBF<S_, 3>(points, distances, new RBFThinPlateSpline<S_>(1.0));
  //rbf = RBF<S_, 3>(points, distances, new RBFBiharmonic<S_>);

  S_ max_error = std::numeric_limits<S_>::min();
  for(unsigned int i = 0; i < points.size(); i++) {
    S_ error = fabs(distances[i] - rbf.eval(points[i]));
    if(error > max_error) max_error = error;
  }

  return max_error;
}

template<class S_, int dim_>
void discretize_rbf(const RBF<S_, dim_> &rbf, ScalarField<S_, dim_> &field) {
  for(int i = 0; i < field.size(); i++) {
    Vector<S_, dim_> x = field.coordinateAt(i); 
    field.at(i) = rbf.eval(x);
  }
}

#endif
