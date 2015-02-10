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
    bool empty();

    T top() const;
    T pop();

    void insert(const std::vector<T> &elements);
    void insert(const T &element);

private:
    unsigned int parent_idx(unsigned int i) const;
    unsigned int left_child_idx(unsigned int i) const;
    unsigned int right_child_idx(unsigned int i) const;
    
    bool comp(const T &a, const T &b);

    void heapify(unsigned int i);

private:
    HeapType m_type;
    std::vector<T> m_data;
};

template <class T>
BinaryHeap<T>::BinaryHeap(HeapType type) : m_type(type) {}

template <class T>
inline
unsigned int BinaryHeap<T>::size() const {
    return m_data.size();
}

template <class T>
inline
T BinaryHeap<T>::top() const {
    return m_data[0];
}

template <class T>
T BinaryHeap<T>::pop() {
    if(m_data.empty())
	throw("Heap underflow");
    
    T max = m_data[0];
    
    m_data[0] = m_data.back();
    m_data.pop_back();
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
    unsigned int e = m_data.size();    

    m_data.push_back(element);
    
    while(e > 0 && comp(m_data[e], m_data[parent_idx(e)])) {
	std::swap(m_data[parent_idx(e)], m_data[e]);
	e = parent_idx(e);
    }
}

template <class T>
inline
bool BinaryHeap<T>::empty() {
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
bool BinaryHeap<T>::comp(const T &a, const T &b) {
    return (m_type == MAX ? a > b : a < b);
}

template <class T>
void BinaryHeap<T>::heapify(unsigned int i) {
    unsigned left = left_child_idx(i), right = right_child_idx(i);
    unsigned int idx = i;
    
    if(left < m_data.size() && comp(m_data[left], m_data[idx]))
	idx = left;
    
    if(right <= m_data.size() && comp(m_data[right], m_data[idx]))
	idx = right;
    
    if(idx != i) {
	std::swap(m_data[i], m_data[idx]);
	heapify(idx);
    }
}

}

#endif
