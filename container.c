#include "llsm.h"
#include <string.h>
#include <stdlib.h>

FP_TYPE* llsm_create_fp(FP_TYPE x) {
  FP_TYPE* ret = malloc(sizeof(FP_TYPE));
  *ret = x;
  return ret;
}

int*     llsm_create_int(int x) {
  int* ret = malloc(sizeof(int));
  *ret = x;
  return ret;
}

FP_TYPE* llsm_create_fparray(int size) {
  void* ret = calloc(sizeof(FP_TYPE) * size + sizeof(int), 1);
  ((int*)ret)[0] = size;
  return (FP_TYPE*)((int*)ret + 1);
}

FP_TYPE* llsm_copy_fp(FP_TYPE* src) {
  return llsm_create_fp(*src);
}

int*     llsm_copy_int(int* src) {
  return llsm_create_int(*src);
}

FP_TYPE* llsm_copy_fparray(FP_TYPE* src) {
  int size = ((int*)src - 1)[0];
  FP_TYPE* ret = llsm_create_fparray(size);
  memcpy(ret, src, size * sizeof(FP_TYPE));
  return ret;
}

void     llsm_delete_fp(FP_TYPE* dst) {
  free(dst);
}

void     llsm_delete_int(int* dst) {
  free(dst);
}

void     llsm_delete_fparray(FP_TYPE* dst) {
  free((int*)dst - 1);
}

int      llsm_fparray_length(FP_TYPE* src) {
  return ((int*)src - 1)[0];
}

llsm_container* llsm_create_container(int nmember) {
  llsm_container* ret = malloc(sizeof(llsm_container));
  ret -> members = calloc(nmember, sizeof(void*));
  ret -> destructors = calloc(nmember, sizeof(llsm_fdestructor));
  ret -> copyctors = calloc(nmember, sizeof(llsm_fcopy));
  ret -> nmember = nmember;
  return ret;
}

llsm_container* llsm_copy_container(llsm_container* src) {
  llsm_container* ret = llsm_create_container(src -> nmember);
  for(int i = 0; i < src -> nmember; i ++) {
    if(src -> copyctors[i] != NULL) {
      ret -> members[i] = src -> copyctors[i](src -> members[i]);
      ret -> destructors[i] = src -> destructors[i];
    } else
      ret -> members[i] = src -> members[i];
    ret -> copyctors[i] = src -> copyctors[i];
  }
  return ret;
}

void llsm_delete_container(llsm_container* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nmember; i ++) {
    if(dst -> destructors[i] != NULL)
      dst -> destructors[i](dst -> members[i]);
  }
  free(dst -> members);
  free(dst -> destructors);
  free(dst -> copyctors);
  free(dst);
}

void* llsm_container_get(llsm_container* src, int index) {
  if(index >= src -> nmember) return NULL;
  return src -> members[index];
}

void llsm_container_attach_(llsm_container* dst, int index, void* ptr,
  llsm_fdestructor dtor, llsm_fcopy copyctor) {
  if(index >= dst -> nmember) { // expand container
    int new_size = index + 1;
    dst -> members = realloc(dst -> members, sizeof(void*) * new_size);
    dst -> destructors = realloc(dst -> destructors,
      sizeof(llsm_fdestructor*) * new_size);
    dst -> copyctors = realloc(dst -> copyctors,
      sizeof(llsm_fcopy*) * new_size);
    for(int i = dst -> nmember; i < new_size; i ++) {
      dst -> members[i] = NULL;
      dst -> destructors[i] = NULL;
      dst -> copyctors[i] = NULL;
    }
    dst -> nmember = new_size;
  }
  llsm_container_remove(dst, index);
  dst -> members[index] = ptr;
  dst -> destructors[index] = dtor;
  dst -> copyctors[index] = copyctor;
}

void llsm_container_remove(llsm_container* dst, int index) {
  if(dst -> members[index] == NULL) return;

  if(dst -> destructors[index] != NULL)
    dst -> destructors[index](dst -> members[index]);
  dst -> members[index] = NULL;
  dst -> destructors[index] = NULL;
  dst -> copyctors[index] = NULL;
}

llsm_chunk* llsm_create_chunk(llsm_container* conf, int init_frames) {
  llsm_chunk* ret = malloc(sizeof(llsm_chunk));
  int* nfrm = llsm_container_get(conf, LLSM_CONF_NFRM);
  int* nchannel = llsm_container_get(conf, LLSM_CONF_NCHANNEL);
  int* npsd = llsm_container_get(conf, LLSM_CONF_NPSD);
  if(nchannel == NULL || npsd == NULL) return NULL;

  ret -> conf = llsm_copy_container(conf);
  if(nfrm != NULL) {
    ret -> frames = calloc(*nfrm, sizeof(llsm_container*));
    if(init_frames)
      for(int i = 0; i < *nfrm; i ++)
        ret -> frames[i] = llsm_create_frame(0, *nchannel, 0, *npsd);
  } else
    ret -> frames = NULL;
  return ret;
}

llsm_chunk* llsm_copy_chunk(llsm_chunk* src) {
  llsm_chunk* ret = llsm_create_chunk(src -> conf, 0);
  int* nfrm = llsm_container_get(src -> conf, LLSM_CONF_NFRM);
  if(nfrm != NULL)
    for(int i = 0; i < *nfrm; i ++)
      ret -> frames[i] = llsm_copy_container(src -> frames[i]);
  return ret;
}

void llsm_delete_chunk(llsm_chunk* dst) {
  if(dst == NULL) return;
  int* nfrm = llsm_container_get(dst -> conf, LLSM_CONF_NFRM);
  if(nfrm != NULL) {
    for(int i = 0; i < *nfrm; i ++)
      llsm_delete_container(dst -> frames[i]);
  }
  llsm_delete_container(dst -> conf);
  free(dst -> frames);
  free(dst);
}
