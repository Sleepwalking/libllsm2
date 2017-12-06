/*
  libllsm2 - Low Level Speech Model (version 2)
  ===
  Copyright (c) 2017 Kanru Hua.

  libllsm2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libllsm2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libllsm. If not, see <http://www.gnu.org/licenses/>.
*/

/** @file */

#ifndef LLSM_BUFFER_H
#define LLSM_BUFFER_H

#include <assert.h>
#include <stdlib.h>

/** @defgroup group_ringbuffer llsm_ringbuffer
 *  @{ */
/** @brief A circular array structure for storing audio samples streaming in
 *    and streaming out in real time. */
typedef struct {
  FP_TYPE* data;  /**< the audio samples */
  int capacity;   /**< the size of the ring buffer */
  int curr;       /**< the current access position */
} llsm_ringbuffer;

/** @brief Create an empty ring buffer with a given size. */
static inline llsm_ringbuffer* llsm_create_ringbuffer(int capacity) {
  assert(capacity > 0);
  llsm_ringbuffer* ret = (llsm_ringbuffer*)malloc(sizeof(llsm_ringbuffer));
  ret -> capacity = capacity;
  ret -> data = (FP_TYPE*)calloc(capacity, sizeof(FP_TYPE));
  ret -> curr = 0;
  return ret;
}

/** @brief Delete and free a ring buffer. */
static inline void llsm_delete_ringbuffer(llsm_ringbuffer* dst) {
  if(dst == NULL) return;
  free(dst -> data);
  free(dst);
}

/** @brief Access the sample at an index relative to the current position;
 *    the index is expected to be negative. */
static inline FP_TYPE llsm_ringbuffer_read(llsm_ringbuffer* src, int idx) {
  assert(idx < 0 && idx >= -src -> capacity);
  return src -> data[(src -> curr + idx + src -> capacity) % src -> capacity];
}

/** @brief Modify the sample at an index relative to the current position;
 *    the index is expected to be negative. */
static inline void llsm_ringbuffer_write(llsm_ringbuffer* dst, int idx, FP_TYPE x) {
  assert(idx < 0 && idx >= -dst -> capacity);
  dst -> data[(dst -> curr + idx + dst -> capacity) % dst -> capacity] = x;
}

/** @brief Append a sample to the ring buffer and push the current position
 *    forward by one sample. */
static inline void llsm_ringbuffer_append(llsm_ringbuffer* dst, FP_TYPE x) {
  dst -> data[dst -> curr] = x;
  dst -> curr = (dst -> curr + 1) % dst -> capacity;
}

/** @brief Push the current position forward by several samples without
 *    modifying the content. */
static inline void llsm_ringbuffer_forward(llsm_ringbuffer* dst, int size) {
  dst -> curr = (dst -> curr + size) % dst -> capacity;
}

/** @brief Load a an array of samples into a destination pointer;
 *    the source starts from a negative index (lag). */
static inline void llsm_ringbuffer_readchunk(llsm_ringbuffer* src,
  int lag, int size, FP_TYPE* dst) {
  assert(size > 0);
  assert(lag + size <= 0);
  assert(lag > -src -> capacity);
  int base = src -> curr + src -> capacity;
  for(int i = 0; i < size; i ++)
    dst[i] = src -> data[(base + lag + i) % src -> capacity];
}

/** @brief Write an array of samples from a source pointer;
 *    the destination starts from a negative index (lag). */
static inline void llsm_ringbuffer_writechunk(llsm_ringbuffer* dst,
  int lag, int size, FP_TYPE* src) {
  assert(size > 0);
  assert(lag + size <= 0);
  assert(lag >= -dst -> capacity);
  int base = dst -> curr + dst -> capacity;
  for(int i = 0; i < size; i ++)
    dst -> data[(base + lag + i) % dst -> capacity] = src[i];
}

/** @brief Add an array of samples to an subset of current samples from a
      source pointer; the subset starts from a negative index (lag). */
static inline void llsm_ringbuffer_addchunk(llsm_ringbuffer* dst,
  int lag, int size, FP_TYPE* src) {
  assert(size > 0);
  assert(lag + size <= 0);
  assert(lag >= -dst -> capacity);
  int base = dst -> curr + dst -> capacity;
  for(int i = 0; i < size; i ++)
    dst -> data[(base + lag + i) % dst -> capacity] += src[i];
}

/** @brief Append an array of samples from a source pointer. This will move
      the current position forward by the size of the data. */
static inline void llsm_ringbuffer_appendchunk(llsm_ringbuffer* dst,
  int size, FP_TYPE* src) {
  assert(size > 0);
  assert(size <= dst -> capacity);
  llsm_ringbuffer_forward(dst, size);
  llsm_ringbuffer_writechunk(dst, -size, size, src);
}

/** @brief Append a chunk of zeros. This will move the current position forward
  *   by the number of zeros. */
static inline void llsm_ringbuffer_appendblank(llsm_ringbuffer* dst, int size) {
  assert(size > 0);
  assert(size <= dst -> capacity);
  llsm_ringbuffer_forward(dst, size);
  int base = dst -> curr + dst -> capacity;
  for(int i = 0; i < size; i ++)
    dst -> data[(base - size + i) % dst -> capacity] = 0;
}
/** @} */

/** @defgroup group_vringbuffer llsm_vringbuffer
 *  @{ */
/** @brief A circular array structure for storing structural objects streaming
 *    int and streaming out in real time. A destructor function has to be
 *    specified (and hence the name virtual ring buffer). */
typedef struct {
  void** data;
  int capacity;
  int curr;
  llsm_fdestructor destructor;
} llsm_vringbuffer;

/** @brief Create an empty virtual ring buffer with a given size. */
static inline llsm_vringbuffer* llsm_create_vringbuffer(int capacity,
  llsm_fdestructor destructor) {
  llsm_vringbuffer* ret = (llsm_vringbuffer*)malloc(sizeof(llsm_vringbuffer));
  ret -> capacity = capacity;
  ret -> curr = 0;
  ret -> data = (void**)malloc(sizeof(void*) * capacity);
  ret -> destructor = destructor;
  for(int i = 0; i < capacity; i ++) ret -> data[i] = NULL;
  return ret;
}

/** @brief Delete and free a virtual ring buffer, including all its stored
 *    objects. */
static inline void llsm_delete_vringbuffer(llsm_vringbuffer* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> capacity; i ++)
    if(dst -> data[i] != NULL)
      dst -> destructor(dst -> data[i]);
  free(dst -> data);
  free(dst);
}

/** @brief Access the object at an index relative to the current position;
 *    the index is expected to be negative. */
static inline void* llsm_vringbuffer_read(llsm_vringbuffer* src, int idx) {
  assert(idx < 0 && idx >= -src -> capacity);
  return src -> data[(src -> curr + idx + src -> capacity) % src -> capacity];
}

/** @brief Replace the object at an index relative to the current position by
 *    a given pointer; no actual deep copy is performed; the index is expected
 *    to be negative. */
static inline void llsm_vringbuffer_write(llsm_vringbuffer* dst, int idx, void* x) {
  assert(idx < 0 && idx >= -dst -> capacity);
  void** currptr = & dst -> data[(dst -> curr + idx + dst -> capacity) % dst -> capacity];
  if(*currptr != NULL) dst -> destructor(*currptr);
  *currptr = x;
}

/** @brief Append an object to the virtual ring buffer and push the current
 *    position forward by one sample; no actual deep copy is performed. */
static inline void llsm_vringbuffer_append(llsm_vringbuffer* dst, void* x) {
  if(dst -> data[dst -> curr] != NULL)
    dst -> destructor(dst -> data[dst -> curr]);
  dst -> data[dst -> curr] = x;
  dst -> curr = (dst -> curr + 1) % dst -> capacity;
}
/** @} */

#endif
