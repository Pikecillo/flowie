#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__

#include <vector>

template <class T>
struct Node;

template <class T, class W>
struct Edge {
  Node<T> *first, *second;
  W weight;
};

template <class T, class W>
struct Node {
  T data;
  
};

template <class T, class W>
class Network {
public:
  std::vector<Node<T, W> > nodes;
  std::vector<Edge<T, W> > edges;

public:
  void addNode(const Node<T, W> &node) {
    nodes.push_back(node);
  }

  void addEdge() {

  }
}

class Digraph {

};

#endif
