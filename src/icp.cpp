/*
 * University of Houston
 * Mario Rincon-Nigro. July 2013.
 */

#include "flowie.hpp"
#include "icp.hpp"
#include "acceleration/kdtree.hpp"
#include "numerical/linear_solver.hpp"

// Registers source to target using an affine transformation
Mat3x4d affine_registration(std::vector<Vector3d> source,
			    std::vector<Vector3d> target) {
  if(source.size() != target.size())
    throw("Source and target have different sizes");

  Matrixd A(source.size(), 4), x(4, 3), b(source.size(), 3);
  
  // Fill matrices to solve linear system
  for(unsigned int i = 0; i < source.size(); i++) {
    Vector3d s = source[i];
    Vector3d t = target[i];
    
    A.set(i, 0, s.x()); A.set(i, 1, s.y());
    A.set(i, 2, s.z()); A.set(i, 3, 1.0);

    b.set(i, 0, t.x()); b.set(i, 1, t.y()); b.set(i, 2, t.z());
  }
  
  // Solve it
  Solver::dgels(A, b, x);
  
  // Copy transformation matrix
  Mat3x4d transf;

  for(int i = 0; i < x.rows(); i++)
    for(int j = 0; j < x.cols(); j++)
      transf[j][i] = x.get(i, j);

  return transf;
}

// Register point cloud source to point cloud target using
// rigid Iterative Closest Point
Mat3x4d icp(std::vector<Vector3d> source,
	    std::vector<Vector3d> target,
	    int max_iter) {
  Mat3x4d rt = Mat3x4d::identity();
  KDTree<3> kdtree(target);
  std::vector<Vector3d> tcorresp;
  int iter = 0;

  while(iter++ < max_iter) {
    // Create correspondences
    tcorresp.clear();
    
    for(unsigned int i = 0; i < source.size(); i++) {
      int ci = kdtree.nearest(source[i]);
      tcorresp.push_back(kdtree.ptss[ci]);
    }

    // Solver for affine transformation
    Mat3x4d curr_rt = affine_registration(source, tcorresp);

    // Transform source point cloud
    for(unsigned int i = 0; i < source.size(); i++)
      source[i] = curr_rt * source[i].homogeneous();

    // Accumulate the transformation
    Mat3x4d temp;
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 4; j++) {
	Vector4d col = rt.col(j).homogeneous();
	if(j != 3) col[3] = 0.0;
	temp[i][j] = curr_rt.row(i).dot(col);
      }

    rt = temp;
  }

  return rt;
}
