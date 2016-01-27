/*
 * University of Houston
 * Mario Rincon-Nigro. May 2013.
 */

#include "geometry.hpp"

template <class T_>
class OctreeNode {
public:
  Vector<3, T_> center;
  T_ halfLength;

public:
  virtual bool isLeaf() const = 0;
};

template <class T_>
class OctreeInnerNode : public OctreeNode<T_> {
public:
  OctreeNode<T_> *children[8];

public:
  virtual bool isLeaf() const { return false; }
};

template <class T_, class P_>
class OctreeLeafNode : public OctreeNode<T_> {
public:
  std::vector<P_> objects;

public:
  virtual bool isLeaf() const { return true; }
};

template <class P, class T_>
class Octree {
public:
  OctreeNode<T_> *root;

public:
  void insert(const P_ &object) {
    
  }

private:
  void insert(OctreeNode<T_> *node, const P_ &object) {
    if(node->isLeaf()) {
      OctreeLeafNode *leaf = static_cast<OctreeLeafNode *>(node);

      
    } else {
      OctreeInnerNode *inner = static_cast<OctreeInnerNode *>(node);

      int l = 0;
      if(l > inner->center)
    }
  }
};

template <class T_, class P_>
OctreeNode<T_> *build_branch(std::vector<P_> &items) {
  
}

template <class T_, class P_>
void build(Octree<T_, P_> &octree, std::vector<P_> &items) {

}
