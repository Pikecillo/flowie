/**
 * Copyright (C) 2013. Mario Rincon-Nigro.
 *
 * This file is a part of Flowie.
 *
 * Flowie is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Flowie is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Flowie.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __BINARY_HEAP_HPP__
#define __BINARY_HEAP_HPP__

#include <algorithm>
#include <functional>
#include <vector>

namespace flowie
{

enum HeapType{ MAX, MIN };

template <class T>
class BinaryHeap
{
public:
    BinaryHeap(HeapType type);

    unsigned int size() const;
    bool empty() const;
    T top() const;
    T pop();
    void insert(const std::vector<T> &elements);
    void insert(const T &element);

private:
    unsigned int parent_idx(unsigned int i) const;
    unsigned int left_child_idx(unsigned int i) const;
    unsigned int right_child_idx(unsigned int i) const;
    
    bool comp(const T &a, const T &b) const;

    void heapify(unsigned int i);

private:
    HeapType m_type;
    std::vector<T> m_items;
};

template <class T>
BinaryHeap<T>::BinaryHeap(HeapType type) : m_type(type) {}

template <class T>
inline
unsigned int BinaryHeap<T>::size() const {
    return m_items.size();
}

template <class T>
inline
T BinaryHeap<T>::top() const {
    return m_items[0];
}

template <class T>
T BinaryHeap<T>::pop() {
    if(m_items.empty())
	throw("Heap underflow");
    
    T max = m_items[0];
    
    m_items[0] = m_items.back();
    m_items.pop_back();
    heapify(0);
    
    return max;
}

template <class T>
void BinaryHeap<T>::insert(const std::vector<T> &elements) {
    for(unsigned int i = 0; i < elements.size(); i++)
	insert(elements[i]);
}

template <class T>
void BinaryHeap<T>::insert(const T &element) {
    unsigned int e = m_items.size();    

    m_items.push_back(element);
    
    while(e > 0 && comp(m_items[e], m_items[parent_idx(e)])) {
	std::swap(m_items[parent_idx(e)], m_items[e]);
	e = parent_idx(e);
    }
}

template <class T>
inline
bool BinaryHeap<T>::empty() const {
    return (size() == 0);
}

template <class T>
inline
unsigned int BinaryHeap<T>::parent_idx(unsigned int i) const {
    return (i - 1) / 2;
}

template <class T>
inline 
unsigned int BinaryHeap<T>::left_child_idx(unsigned int i) const {
    return 2 * i + 1;
}

template <class T>
inline
unsigned int BinaryHeap<T>::right_child_idx(unsigned int i) const {
    return 2 * i + 2;
}

template <class T>
bool BinaryHeap<T>::comp(const T &a, const T &b) const {
    return (m_type == MAX ? a > b : a < b);
}

template <class T>
void BinaryHeap<T>::heapify(unsigned int i) {
    unsigned left = left_child_idx(i), right = right_child_idx(i);
    unsigned int idx = i;
    
    if(left < m_items.size() && comp(m_items[left], m_items[idx]))
	idx = left;
    
    if(right <= m_items.size() && comp(m_items[right], m_items[idx]))
	idx = right;
    
    if(idx != i) {
	std::swap(m_items[i], m_items[idx]);
	heapify(idx);
    }
}

}

#endif
