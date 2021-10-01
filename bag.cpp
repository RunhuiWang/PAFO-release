#ifndef BAG_CPP
#define BAG_CPP

#include <stdlib.h>
#include <string.h>
#include "bag.h"

#define typeof(x) __typeof__(x)
// static inline uint
// MAX(uint a, uint b) { return a ^ ( (a ^ b) & -(a < b) ); }
#define MAX(a, b) \
  ({ typeof(a) _a = (a); \
    typeof(b) _b = (b); \
    _a ^ ( (_a ^ _b) & -(_a < _b) ); })

//////////////////////////////////
///                            ///
/// Pennant method definitions ///
///                            ///
//////////////////////////////////
template <typename T>
Pennant<T>::Pennant()
{
  this->els = new T[BLK_SIZE];
  this->l = NULL;
  this->r = NULL;
}

// els_array must have size BLK_SIZE
template <typename T>
Pennant<T>::Pennant(T els_array[])
{
  this->els = els_array;
  this->l = NULL;
  this->r = NULL;
}

template <typename T>
Pennant<T>::~Pennant()
{
  delete[] els;
}

template <typename T>
inline const T*
Pennant<T>::getElements()
{
  return this->els;
}

template <typename T>
inline Pennant<T>*
Pennant<T>::getLeft()
{
  return l;
}

template <typename T>
inline Pennant<T>*
Pennant<T>::getRight()
{
  return r;
}

/*
 * This method assumes that the pennant rooted at <that> has the
 * same number of nodes as the pennant rooted at <this>.
 */
template <typename T>
inline Pennant<T>*
Pennant<T>::combine(Pennant<T>* that)
{
  that->r = this->l;
  this->l = that;

  return this;
}

/*
 * This method performs the inverse operation of
 * Pennant<T>::combine and places half of the split here and the
 * the other half at the returned pointer
 */
template <typename T>
inline Pennant<T>*
Pennant<T>::split()
{
  Pennant<T>* that;

  that = this->l;
  this->l = that->r;
  that->r = NULL;

  return that;
}

//////////////////////////////
///                        ///
/// Bag method definitions ///
///                        ///
//////////////////////////////
template <typename T>
Bag<T>::Bag() : fill(0), size(0)
{
  this->bag = new Pennant<T>*[BAG_SIZE];
#if FILLING_ARRAY
  this->filling = new T[BLK_SIZE];
#else
  this->filling = new Pennant<T>();
#endif  // FILLING_ARRAY

}

/*
 * Copy Constructor. Performs a shallow bag copy.
 * Useful for verifying Bag implementation correctness, since
 * the split method destroys the Bag.
 */
template <typename T>
Bag<T>::Bag(Bag<T> *that) : fill(that->fill), size(that->size)
{
  this->bag = new Pennant<T>*[BAG_SIZE];
  for (uint i = 0; i < BAG_SIZE; i++)
    this->bag[i] = that->bag[i];

  this->filling = that->filling;
}

template <typename T>
Bag<T>::~Bag()
{
  if (this->bag != NULL) {
    for (int i = 0; i < this->fill; i++) {
      if (this->bag[i] != NULL) {
	delete this->bag[i];
	this->bag[i] = NULL;
      }
    }

    delete[] this->bag;
  }
#if FILLING_ARRAY
  delete[] this->filling;
#else
  delete this->filling;
#endif  // FILLING_ARRAY
}

template <typename T>
inline uint
Bag<T>::numElements() const
{
  uint count = this->size;
  uint k = 1;
  for (uint i = 0; i < this->fill; i++) {
    if (this->bag[i] != NULL)
      count += k*BLK_SIZE;
    k = k * 2;
  }

  return count;
}

// helper routine to perform bag-insert with filled filling array
template <typename T>
inline void
Bag<T>::insert_h()
{
#if FILLING_ARRAY
  Pennant<T> *c = new Pennant<T>(this->filling);
  this->filling = new T[BLK_SIZE];
#else
  Pennant<T> *c = this->filling;
  this->filling = new Pennant<T>();
#endif  // FILLING_ARRAY
  //this->filling->els[0] = el;
  //this->size = 1;
  this->size = 0;

  uint i = 0;
  do {

    if (i < this->fill && this->bag[i] != NULL) {

      c = this->bag[i]->combine(c);
      this->bag[i] = NULL;

    } else {

      this->bag[i] = c;
      //this->fill = this->fill > i+1 ? this->fill : i+1;
      this->fill = MAX(this->fill, i+1);
      return;

    }
    // switch(i < this->fill && this->bag[i] != NULL) {
    // case false:
    //   this->bag[i] = c;
    //   this->fill = this->fill > i+1 ? this->fill : i+1;
    //   return;
    // case true:
    //   c = this->bag[i]->combine(c);
    //   this->bag[i] = NULL;
    //   break;
    // default: break;
    // }

    ++i;

  } while (i < BAG_SIZE);

  this->fill = BAG_SIZE;
}

// helper routine to perform bag-insert with filled filling array
template <typename T>
inline void
Bag<T>::insert_fblk(T fblk[])
{
  Pennant<T> *c = new Pennant<T>(fblk);

  uint i = 0;
  do {

    if (i < this->fill && this->bag[i] != NULL) {

      c = this->bag[i]->combine(c);
      this->bag[i] = NULL;

    } else {

      this->bag[i] = c;
      //this->fill = this->fill > i+1 ? this->fill : i+1;
      this->fill = MAX(this->fill, i+1);
      return;

    }
    // switch(i < this->fill && this->bag[i] != NULL) {
    // case false:
    //   this->bag[i] = c;
    //   this->fill = this->fill > i+1 ? this->fill : i+1;
    //   return;
    // case true:
    //   c = this->bag[i]->combine(c);
    //   this->bag[i] = NULL;
    //   break;
    // default: break;
    // }

    ++i;

  } while (i < BAG_SIZE);

  this->fill = BAG_SIZE;
}

template <typename T>
inline void
Bag<T>::insert_blk(T blk[], uint size)
{
  int i;

  // Deal with the partially-filled Pennants
  if (this->size < size) {
    // Copy contents of this->filling into blk
    i = this->size - (BLK_SIZE - size);

    if (i >= 0) {
      // Contents of this->filling fill blk
#if FILLING_ARRAY
      memcpy(blk + size,
             this->filling + i,
             (BLK_SIZE - size) * sizeof(T));
#else
      memcpy(blk + size,
             this->filling->els + i,
             (BLK_SIZE - size) * sizeof(T));
#endif  // FILLING_ARRAY

      //carry = blk;
      insert_fblk(blk);

      this->size = i;
    } else {
      // Contents of this->filling do not fill blk
      //if (this->size > 0) {
#if FILLING_ARRAY
      memcpy(blk + size,
             this->filling,
             this->size * sizeof(T));
#else
      memcpy(blk + size,
             this->filling->els,
             this->size * sizeof(T));
#endif  // FILLING_ARRAY

      delete[] this->filling;

#if FILLING_ARRAY
      this->filling = blk;
#else
      this->filling->els = blk;
#endif  // FILLING_ARRAY
      this->size += size;
    }

  } else {
    // Copy contents of blk into this->filling
    T* carry;
    i = size - (BLK_SIZE - this->size);

    if (i >= 0) {
      // Contents of blk fill this->filling
#if FILLING_ARRAY
      memcpy(this->filling + this->size,
             blk + i,
             (BLK_SIZE - this->size) * sizeof(T));

      carry = this->filling;
      this->filling = blk;

#else
      memcpy(this->filling->els + this->size,
             blk + i,
             (BLK_SIZE - this->size) * sizeof(T));

      carry = this->filling->els;
      this->filling->els = blk;

#endif  // FILLING_ARRAY

      this->size = i;
      insert_fblk(carry);

    } else {
      // Contents of blk do not fill this->filling
      //if (that->size > 0) {
#if FILLING_ARRAY
	memcpy(this->filling + this->size,
	       blk,
	       size * sizeof(T));
#else
	memcpy(this->filling->els + this->size,
	       blk,
	       size * sizeof(T));
#endif  // FILLING_ARRAY
        //}
      this->size += size;
    }
  }
  
}

template <typename T>
inline void
Bag<T>::insert(T el)
{
  //assert(this->size < BLK_SIZE);
#if FILLING_ARRAY
  this->filling[this->size++] = el;
#else
  this->filling->els[this->size++] = el;
#endif  // FILLING_ARRAY

  if (this->size < BLK_SIZE) {
    return;
  }

#if FILLING_ARRAY
  Pennant<T> *c = new Pennant<T>(this->filling);
  this->filling = new T[BLK_SIZE];
#else
  Pennant<T> *c = this->filling;
  this->filling = new Pennant<T>();
#endif  // FILLING_ARRAY
  //this->filling->els[0] = el;
  //this->size = 1;
  this->size = 0;

  int i = 0;
  do {

    if (i < this->fill && this->bag[i] != NULL) {

      c = this->bag[i]->combine(c);
      this->bag[i] = NULL;

    } else {

      this->bag[i] = c;
      //this->fill = this->fill > i+1 ? this->fill : i+1;
      this->fill = MAX(this->fill, i+1);
      return;

    }

    ++i;

  } while (i < BAG_SIZE);

  this->fill = BAG_SIZE;
}

template <typename T>
void
Bag<T>::merge(Bag<T>* that)
{

  Pennant<T> *c = NULL;
#if FILLING_ARRAY
  T* carry = NULL;
#endif  // FILLING_ARRAY
  char x;
  int i;

  // Deal with the partially-filled Pennants
  if (this->size < that->size) {
    i = this->size - (BLK_SIZE - that->size);

    if (i >= 0) {
#if FILLING_ARRAY
      memcpy(that->filling + that->size,
             this->filling + i,
             (BLK_SIZE - that->size) * sizeof(T));

      carry = that->filling;
#else
      memcpy(that->filling->els + that->size,
             this->filling->els + i,
             (BLK_SIZE - that->size) * sizeof(T));

      c = that->filling;
#endif  // FILLING_ARRAY
      this->size = i;
    } else {
      //if (this->size > 0) {
#if FILLING_ARRAY
	memcpy(that->filling + that->size,
	       this->filling,
	       this->size * sizeof(T));
#else
        memcpy(that->filling->els + that->size,
	       this->filling->els,
	       this->size * sizeof(T));
#endif  // FILLING_ARRAY
        //}
#if FILLING_ARRAY
      delete[] this->filling;
#else
      delete this->filling;
#endif  // FILLING_ARRAY

      this->filling = that->filling;
      this->size += that->size;
    }

  } else {
    i = that->size - (BLK_SIZE - this->size);

    if (i >= 0) {

#if FILLING_ARRAY
      memcpy(this->filling + this->size,
             that->filling + i,
             (BLK_SIZE - this->size) * sizeof(T));

      carry = this->filling;
#else
      memcpy(this->filling->els + this->size,
             that->filling->els + i,
             (BLK_SIZE - this->size) * sizeof(T));

      c = this->filling;
#endif  // FILLING_ARRAY

      this->filling = that->filling;
      this->size = i;
    } else {

      //if (that->size > 0) {
#if FILLING_ARRAY
	memcpy(this->filling + this->size,
	       that->filling,
	       that->size * sizeof(T));
#else
	memcpy(this->filling->els + this->size,
	       that->filling->els,
	       that->size * sizeof(T));
#endif  // FILLING_ARRAY
        //}
      this->size += that->size;
    }
  }

  that->filling = NULL;

  // Update this->fill (assuming no final carry)
  int min, max;
  if (this->fill > that->fill) {
    min = that->fill;
    max = this->fill;
  } else {
    min = this->fill;
    max = that->fill;
  }

#if FILLING_ARRAY
  if (carry != NULL)
    c = new Pennant<T>(carry);
#endif  // FILLING_ARRAY

  // Merge
  for (i = 0; i < min; ++i) {

    x =
      (this->bag[i] != NULL) << 2 |
      (i < that->fill && that->bag[i] != NULL) << 1 |
      (c != NULL);

    switch(x) {
    case 0x0:
      this->bag[i] = NULL;
      c = NULL;
      that->bag[i] = NULL;
      break;
    case 0x1:
      this->bag[i] = c;
      c = NULL;
      that->bag[i] = NULL;
      break;
    case 0x2:
      this->bag[i] = that->bag[i];
      that->bag[i] = NULL;
      c = NULL;
      break;
    case 0x3:
      c = that->bag[i]->combine(c);
      that->bag[i] = NULL;
      this->bag[i] = NULL;
      break;
    case 0x4:
      /* this->bag[i] = this->bag[i]; */
      c = NULL;
      that->bag[i] = NULL;
      break;
    case 0x5:
      c = this->bag[i]->combine(c);
      this->bag[i] = NULL;
      that->bag[i] = NULL;
      break;
    case 0x6:
      c = this->bag[i]->combine(that->bag[i]);
      that->bag[i] = NULL;
      this->bag[i] = NULL;
      break;
    case 0x7:
      /* this->bag[i] = this->bag[i]; */
      c = that->bag[i]->combine(c);
      that->bag[i] = NULL;
      break;
    default: break;
    }
  }

  that->fill = 0;

  if (this->fill == max) {
    if (c == NULL)
      return;

    do {
      if (i < max && this->bag[i] != NULL) {

	c = this->bag[i]->combine(c);
	this->bag[i] = NULL;

      } else {

	this->bag[i] = c;
	//this->fill = max > i+1 ? max : i+1;
        this->fill = MAX(max, i+1);
	return;

      }
      // switch(i < max && this->bag[i] != NULL) {
      // case false:
      //   this->bag[i] = c;
      //   this->fill = max > i+1 ? max : i+1;
      //   return;
      // case true:
      //   c = this->bag[i]->combine(c);
      //   this->bag[i] = NULL;
      //   break;
      // default: break;
      // }
      
      ++i;
      
    } while (i < BAG_SIZE);
  } else { // that->fill == max
    int j;
    if (c == NULL) {
      this->fill = max;
      for (j = i; j < this->fill; ++j) {
	this->bag[j] = that->bag[j];
	that->bag[j] = NULL;
      }
      return;
    }

    do {
      if (i < max && that->bag[i] != NULL) {

	c = that->bag[i]->combine(c);
	this->bag[i] = NULL;
	that->bag[i] = NULL;

      } else {

	this->bag[i] = c;
	//this->fill = max > i+1 ? max : i+1;
        this->fill = MAX(max, i+1);

	for (j = i+1; j < this->fill; ++j) {
	  this->bag[j] = that->bag[j];
	  that->bag[j] = NULL;
	}
	return;

      }
      // switch(i < max && that->bag[i] != NULL) {
      // case false:
      //   this->bag[i] = c;
      //   this->fill = max > i+1 ? max : i+1;
      //   for (j = i+1; j < this->fill; ++j) {
      //     this->bag[j] = that->bag[j];
      //     that->bag[j] = NULL;
      //   }
      //   return;
      // case true:
      //   c = that->bag[i]->combine(c);
      //   this->bag[i] = NULL;
      //   that->bag[i] = NULL;
      //   break;
      // default: break;
      // }
      
      ++i;
      
    } while (i < BAG_SIZE);    
  }

  this->fill = BAG_SIZE;
}

template <typename T>
inline bool
Bag<T>::split(Pennant<T>* *p)
{
  if (this->fill == 0) {
    *p = NULL;
    return false;
  }

  this->fill--;
  *p = this->bag[this->fill];
  this->bag[this->fill] = NULL;

  for ( /* this->fill */; this->fill > 0; this->fill--) {
    if (this->bag[this->fill-1] != NULL)
      break;
  }
  return true;
}

template <typename T>
int
Bag<T>::split(Pennant<T>* *p, int pos)
{
  if (pos < 0 || this->fill <= pos) {
    *p = NULL;
    return this->fill - 1;
  }

  *p = this->bag[pos];

  for (int i = pos-1; i >= 0; i--) {
    if (this->bag[i] != NULL)
      return i;
  }

  return -1;
}

template <typename T>
inline uint
Bag<T>::getFill() const
{
  return this->fill;
}

template <typename T>
inline bool
Bag<T>::isEmpty() const
{
  return this->fill == 0 && this->size == 0;
}

template <typename T>
inline Pennant<T>*
Bag<T>::getFirst() const
{
  return this->bag[0];
}

template <typename T>
inline T*
Bag<T>::getFilling() const
{
#if FILLING_ARRAY
  return this->filling;
#else
  return this->filling->els;
#endif  // FILLING_ARRAY
}

template <typename T>
inline uint
Bag<T>::getFillingSize() const
{
  return this->size;
}

template <typename T>
inline void
Bag<T>::clear()
{
  this->fill = 0;
  this->size = 0;
}

#endif  // BAG_CPP
