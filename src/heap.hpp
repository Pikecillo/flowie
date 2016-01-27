/*
 * University of Houston
 * Mario Rincon-Nigro. July 2013.
 */

#include <algorithm>
#include <vector>

/*
 * Maximum binary heap
 */
template <class T>
class BinaryHeap {
public:
  BinaryHeap() {}

  void build(const std::vector<T> &elements) {
    heap = elements;

    for(unsigned int i = heap.size() / 2; i >= 0; i--)
      heapify(i);
  }

  unsigned int size() const { return heap.size(); }
  T maximum() const { return heap[0]; }
  
  T extract_max() {
    if(heap.empty()) throw("Heap underflow");

    T max = heap[0];

    heap[0] = heap.back();
    heap.pop_back();
    heapify(0);

    return max;
  }

  void insert(const T &element) {
    unsigned int e = heap.size();

    heap.push_back(element);

    while(parent(e) >= 0 && heap[parent(e)] < heap[e]) {
      std::swap(heap[parent(e)], heap[e]);
      e = parent(e);
    }
  }

private:
  unsigned int parent(unsigned int i) const { return (i - 1) / 2; }
  unsigned int left(unsigned int i) const { return 2 * i + 1; }
  unsigned int right(unsigned int i) const { return 2 * i + 2; }

  void heapify(unsigned int i) {
    unsigned l = left(i), r = right(i);
    unsigned int largest;

    if(l < heap.size() && heap[i] < heap[l])
      largest = l;
    else
      largest = i;

    if(r <= heap.size() && heap[largest] < heap[r])
      largest = r;

    if(largest != i) {
      std::swap(heap[i], heap[largest]);
      heapify(largest);
    }
  }

  std::vector<T> heap;
};
