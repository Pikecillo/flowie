/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __FLOWIE_HPP__
#define __FLOWIE_HPP__

#include "field_types.hpp"
#include "geometry.hpp"
#include "heap.hpp"
#include "histogram.hpp"
#include "icp.hpp"
#include "matrix_types.hpp"
#include "primitives.hpp"
#include "reconstruction.hpp"
#include "vector_types.hpp"

// Acceleration module
#include <acceleration/aabb.hpp>
#include <acceleration/bvh.hpp>
#include <acceleration/kdtree.hpp>

// Implicit surface module
#include <implicit/blobs.hpp>
#include <implicit/distance_field.hpp>
#include <implicit/dynamic_field.hpp>
#include <implicit/implicit_surface.hpp>
#include <implicit/velocity.hpp>

// Mesh module
#include <mesh/triangle_mesh.hpp>
#include <mesh/triangle_soup.hpp>

// Numerical module
#include <numerical/derivative.hpp>
#include <numerical/interpolant.hpp>
#include <numerical/linear_solver.hpp>
#include <numerical/lmarquardt.hpp>
#include <numerical/mvopt.hpp>
#include <numerical/pca.hpp>
#include <numerical/qnewton.hpp>
#include <numerical/rbf.hpp>

// Utils module
#include <utils/general.hpp>
#include <utils/misc.hpp>
#include <utils/reference.hpp>
#include <utils/timer.hpp>

// Visualization module
#include <visualization/color_scheme.hpp>
#include <visualization/graphics.hpp>
#include <visualization/isolevel.hpp>

#endif
