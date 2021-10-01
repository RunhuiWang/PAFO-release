// -*- C++ -*-
// Copyright (c) 2010, Tao B. Schardl
/*
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#ifndef BAG_H
#define BAG_H

#include <stdlib.h>
#include <assert.h>
#include <cilk/cilk.h>
#include <cilk/reducer.h>

#include <cilk/cilk_api.h>

#define FILLING_ARRAY true

const int BAG_SIZE = 64;
//const int BLK_SIZE = 128;
const int BLK_SIZE = 128;

template <typename T> class Bag;
template <typename T> class Bag_reducer;

template <typename T>
class Pennant
{
private:
  T* els;
  Pennant<T> *l, *r;
  
public:
  Pennant();
  Pennant(T*);
  ~Pennant();
  
  inline const T* getElements();
  inline Pennant<T>* getLeft();
  inline Pennant<T>* getRight();
  
  inline void clearChildren() {
    l = NULL;
    r = NULL;
  }

  inline Pennant<T>* combine(Pennant<T>*);
  inline Pennant<T>* split();
  
  friend class Bag<T>;
  friend class Bag_reducer<T>;
};

template <typename T>
class Bag
{
//private:
public: // HACKS

  // fill points to the bag entry one position beyond the MSP
  int fill;
  Pennant<T>* *bag;

#if FILLING_ARRAY
  T* filling;
#else
  Pennant<T>* filling;
#endif // FILLING_ARRAY

  // size points to the filling array entry on position
  // beyond last valid element.
  uint size;

  inline void insert_h();
  inline void insert_fblk(T* fblk);
  inline void insert_blk(T* blk, uint size);

public:
  Bag();
  Bag(Bag<T>*);
  
  ~Bag();
  
  inline void insert(T);
  void merge(Bag<T>*);
  inline bool split(Pennant<T>**);
  int split(Pennant<T>**, int);
  
  inline uint numElements() const;
  inline uint getFill() const;
  inline bool isEmpty() const;
  inline Pennant<T>* getFirst() const;
  inline T* getFilling() const;
  inline uint getFillingSize() const;
  
  void clear();
  
  friend class Bag_reducer<T>;
  friend class cilk::monoid_base< Bag<T> >;
};

template <typename T>
class Bag_reducer
{
public:
  struct Monoid: cilk::monoid_base< Bag<T> >
  {
    static void reduce (Bag<T> *left, Bag<T> *right) {
      left->merge(right);
    }
  };
  
private:
  cilk::reducer<Monoid> imp_;
    
public:
  Bag_reducer();
    
  void insert(T);
  void merge(Bag_reducer<T>*);
  inline bool split(Pennant<T>**);
  int split(Pennant<T>**, int);
  
  inline Bag<T> &get_reference();
  
  inline uint numElements() const;
  inline uint getFill() const;
  inline bool isEmpty() const;
  inline Pennant<T>* getFirst() const;
  inline T* getFilling() const;
  inline uint getFillingSize() const;
  
  inline void clear();
};

template <typename T>
Bag_reducer<T>::Bag_reducer() : imp_() { }

template <typename T>
void
Bag_reducer<T>::insert(T el)
{
  imp_.view().insert(el);
}

template <typename T>
void
Bag_reducer<T>::merge(Bag_reducer<T>* that)
{
  this->imp_.view().merge(that->imp_.view());
}

template <typename T>
inline bool
Bag_reducer<T>::split(Pennant<T>* *p)
{
  return imp_.view().split(p);
}

template <typename T>
int
Bag_reducer<T>::split(Pennant<T>* *p, int pos)
{
  return imp_.view().split(p, pos);
}

template <typename T>
inline Bag<T>&
Bag_reducer<T>::get_reference()
{
  return imp_.view();
}

template <typename T>
inline uint
Bag_reducer<T>::numElements() const
{
  return imp_.view().numElements();
}

template <typename T>
inline uint
Bag_reducer<T>::getFill() const
{
  return imp_.view().getFill();
}

template <typename T>
inline bool
Bag_reducer<T>::isEmpty() const
{
  return imp_.view().isEmpty();
}

template <typename T>
inline Pennant<T>*
Bag_reducer<T>::getFirst() const
{
  return imp_.view().getFirst();
}

template <typename T>
inline T*
Bag_reducer<T>::getFilling() const
{
  return imp_.view().getFilling();
}

template <typename T>
inline uint
Bag_reducer<T>::getFillingSize() const
{
  return imp_.view().getFillingSize();
}

template <typename T>
inline void
Bag_reducer<T>::clear()
{
  imp_.view().clear();
}

#include "bag.cpp"

#endif
