/*
 * University of Houston
 * Mario Rincon-Nigro. April 2013.
 */

#ifndef __BVH_HPP__
#define __BVH_HPP__

#include <algorithm>
#include <iostream>
#include <limits>

#include "aabb.hpp"
#include "geometry.hpp"
#include "primitives.hpp"
#include "vector_types.hpp"

template <class T_=float>
class BVHNode {
public:
  virtual bool isLeaf() const = 0;

  const AABB<T_, 3> &getAABB() const { return aabb; }
  void setAABB(const AABB<T_, 3> &bb) { aabb = bb; }

private:
  AABB<T_, 3> aabb;
};

template <class T_=float>
class BVHInnerNode : public BVHNode<T_> {
public:
  BVHInnerNode(BVHNode<T_> *l, BVHNode<T_> *r) :
    left(l), right(r) {}

  BVHInnerNode(const BVHInnerNode<T_> &other) {
    (*this) = other;
  }
  
  BVHInnerNode<T_> &operator=(const BVHInnerNode<T_> &other) {
    left = other.left;
    right = other.right;  
    return (*this);
  }

  BVHNode<T_> *getLeft() const { return left; }
  BVHNode<T_> *getRight() const { return right; }

  virtual bool isLeaf() const { return false; }

private:
  BVHNode<T_> *left;
  BVHNode<T_> *right;
};

template <class T_=float>
class BVHLeafNode : public BVHNode<T_> {
public:
  BVHLeafNode(unsigned int s, unsigned int e)
    : start(s), end(e) {}
  
  BVHLeafNode(const BVHLeafNode<T_> &other) {
    (*this) = other;
  }

  BVHLeafNode<T_> &operator=(const BVHLeafNode<T_> &other) {
    start = other.start;
    end = other.end;
    return (*this);
  }

  virtual bool isLeaf() const { return true; }
  int getStart() const { return start; }
  int getEnd() const { return end; }

private:
  unsigned int start;
  unsigned int end;
};

template <class Primitive_, typename T_=float>
class BVH {
public:
  BVH() : root(0) {}

  BVH(const BVH<Primitive_, T_> &other) {
    (*this) = other;
  }

  ~BVH() {
    deallocate(root);
  }

  void operator=(const BVH<Primitive_, T_> &other) {
    primitives = other.primitives;
    root = other.root;
  }

  std::vector<Primitive_> getPrimitives() const { return primitives; }
  void setPrimitives(std::vector<Primitive_> &p) { primitives = p; }

  BVHNode<T_> *getRoot() const { return root; }
  void setRoot(BVHNode<T_> *r) { root = r; }

  Vector<T_, 3> nearestPoint(BVHNode<T_> *node, const Vector<T_, 3> &p) {
    T_ dist1, dist2;
    Vector<T_, 3> nearest1(std::numeric_limits<T_>::max());
    Vector<T_, 3> nearest2(std::numeric_limits<T_>::max());
    BVHNode<T_> *node1, *node2;

     // If the node is a leaf check distance to all primitives
    if(node->isLeaf()) {
      Vector<T_, 3> nearest(std::numeric_limits<T_>::max());
      BVHLeafNode<T_> *leaf = dynamic_cast<BVHLeafNode<T_> *>(node);

      // Go through each primitive
      for(int i = leaf->getStart(); i < leaf->getEnd(); i++) {
	Vector<T_, 3> current = primitives[i].nearestPoint(p);
	if((current - p).length() < (nearest - p).length())
	  nearest = current;
      }

      return nearest;
    }
    
    BVHInnerNode<T_> *inner = dynamic_cast<BVHInnerNode<T_> *>(node);
    node1 = inner->getLeft();
    node2 = inner->getRight();
    dist1 = min_distance(node1->getAABB(), p);
    dist2 = min_distance(node2->getAABB(), p);
  
    if(dist1 > dist2) {
      std::swap(node1, node2);
      std::swap(dist1, dist2);
    }

    // Try nearest child
    nearest1 = nearestPoint(node1, p);

    // Try the other children if necessary
    if((p - nearest1).length() > dist2)
      nearest2 = nearestPoint(node2, p);

    if((p - nearest1).length() < (p - nearest2).length())
      return nearest1;
    return nearest2;
  }

  Vector<T_, 3> nearestPoint(const Vector<T_, 3> &p) {
    return nearestPoint(root, p);
  }

private:
  void deallocate(BVHNode<T_> *node) {
    if(!node->isLeaf()) {
      deallocate(dynamic_cast<BVHInnerNode<T_> *>(node)->getLeft());
      deallocate(dynamic_cast<BVHInnerNode<T_> *>(node)->getRight());
    }
    
    delete node;
    node = 0;
  }

private:
  std::vector<Primitive_> primitives;
  BVHNode<T_> *root;
};

template <class Primitive_, class T_>
unsigned int split(std::vector<Primitive_> &primitives,
                   unsigned int start, unsigned int end,
                   T_ pivot, unsigned int axis) {
  unsigned int split_index = start;
 
  // Every primitive in [0,split_point) has aabb centroid less than pivot
  // and everything in [split_point,i) has aabb centoird greater than pivot
  for(unsigned int i = start; i < end; i++) {
    AABB<T_, 3> aabb;
    Primitive_ p = primitives[i];

    aabb.join(p.getAABB());
 
    if(aabb.centroid()[axis] < pivot) {
      primitives[i] = primitives[split_index];
      primitives[split_index] = p;
      
      split_index++;
    }
  }
  
  // This ensures the interval is split even in the
  // case all primitives were to one side of the pivot
  if(split_index == start || split_index == end)
    split_index = (start + end) / 2;
  
  return split_index;
}

template <class Primitive_, class T_>
BVHNode<T_> *build_branch(std::vector<Primitive_> &primitives,
			  unsigned int start, unsigned int end,
			  unsigned int axis) {
  AABB<T_, 3> aabb;
  BVHNode<T_> *node = 0;
  const unsigned int max_primitives = 8;

  assert(start <= end);

  // Compute the AABB for this node
  for(unsigned int i = start; i < end; i++)
    aabb.join(primitives[i].getAABB());

  // Too few primitives, make a leaf
  if(end - start <= max_primitives) {
    node = new BVHLeafNode<T_>(start, end);
  } 
  // Otherwise make another inner node
  else {
    T_ pivot = aabb.centroid()[axis];
    int split_index = split<Primitive_, T_>(primitives, start, end, pivot, axis);
    
    axis = (axis + 1) % 3;
    
    BVHNode<T_> *left = build_branch<Primitive_, T_>(primitives, start,
						     split_index, axis);
    BVHNode<T_> *right = build_branch<Primitive_, T_>(primitives, split_index,
						      end, axis);

    node = new BVHInnerNode<T_>(left, right);
  }
  
   // Set the AABB for the current node
  node->setAABB(aabb);

  return node;
}

template <class Primitive_, class T_>
void centroid_build(BVH<Primitive_, T_> &bvh,
		    const std::vector<Primitive_> &primitives) {
  std::vector<Primitive_> p = primitives;

  bvh.setRoot(build_branch<Primitive_, T_>(p, 0, p.size(), 0));
  bvh.setPrimitives(p);
}

#endif
