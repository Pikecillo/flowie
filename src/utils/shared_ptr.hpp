/**
 * Copyright (C) 2014. Mario Rincon-Nigro.
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

#ifndef __SHARED_PTR_HPP__
#define __SHARED_PTR_HPP__

#include <string.h>

namespace flowie
{

template <class T>
class SharedMem
{
public:
    explicit SharedMem(T* ptr);
    ~SharedMem();
    
    void operator=(const SharedMem<T>& other);
    
    long int& use_count();
    T* get() const;
    
private:
    long int m_use_count;
    T* m_ptr;
};

template <class T>
SharedMem<T>::SharedMem(T* ptr) : m_ptr(ptr) {
}

template <class T>
SharedMem<T>::~SharedMem() {
    delete m_ptr;
}

template <class T>
long int& SharedMem<T>::use_count() {
    return m_use_count;
}

template <class T>
T* SharedMem<T>::get() const {
    return m_ptr;
}

template <class T>
class SharedPtr
{
public:
    SharedPtr();
    explicit SharedPtr(T* ptr);
    SharedPtr(const SharedPtr<T> &other);
    ~SharedPtr();
    
    void operator=(const SharedPtr<T> &other);
    
    T& operator*() const;
    T* operator->() const;
    
    T* get() const;
    long int use_count() const;
    bool unique() const;
    void reset(T* ptr);
    void swap(SharedPtr<T> &other);
    void clone(const SharedPtr<T> &other);
    
private:
    void release();
    
private:
    SharedMem<T>* m_shared_mem;
};

template <class T>
SharedPtr<T>::SharedPtr() : m_shared_mem(0) {}

template <class T>
SharedPtr<T>::SharedPtr(T* ptr) : m_shared_mem(0) {
    reset(ptr);
}

template <class T>
SharedPtr<T>::SharedPtr(const SharedPtr<T> &other) {
    (*this) = other;
}

template <class T>
SharedPtr<T>::~SharedPtr() {
    release();
}

template <class T>
void SharedPtr<T>::release() {
    if (m_shared_mem) {
	if (unique()) {
	    delete m_shared_mem;
	}
	else {
	    m_shared_mem->use_count()--;
	}
	
	m_shared_mem = 0;
    }
}

template <class T>
void SharedPtr<T>::operator=(const SharedPtr<T> &other) {
    if (m_shared_mem)
	release();
    
    if (m_shared_mem = other.m_shared_mem) {
	m_shared_mem->use_count()++;
    }
}

template <class T>
T& SharedPtr<T>::operator*() const {
    return *(m_shared_mem->get());
}

template <class T>
T* SharedPtr<T>::operator->() const {
    return m_shared_mem->get();
}

template <class T>
T* SharedPtr<T>::get() const {
    if (m_shared_mem)
	return m_shared_mem->get();
    
    return 0;
}

template <class T>
long int SharedPtr<T>::use_count() const {
    if (m_shared_mem)
	return m_shared_mem->use_count();
    
    return 0;
}

template <class T>
bool SharedPtr<T>::unique() const {
    return (m_shared_mem->use_count() < 2);
}

template <class T>
void SharedPtr<T>::reset(T* ptr) {
    if (m_shared_mem)
	release();
    
    if (ptr)
	m_shared_mem = new SharedMem<T>(ptr);
}

template <class T>
void SharedPtr<T>::swap(SharedPtr<T> &other) {
    SharedMem<T>* t = m_shared_mem;
    m_shared_mem = other.m_shared_mem;
    other.m_shared_mem = m_shared_mem;
}

template <class T>
void SharedPtr<T>::clone(const SharedPtr<T> &other) {
    T* ptr = new T();
    memcpy(ptr, other.get(), sizeof(T));
    reset(ptr);
}

}

#endif
