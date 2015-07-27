/// \file
/// Wrappers for memory allocation.

#ifndef _MEMUTILS_H_
#define _MEMUTILS_H_

#include <stdlib.h>
#include <stdio.h>

#define freeMe(s,element) {if(s->element) comdFree(s->element);  s->element = NULL;}

static inline void* comdMalloc(size_t iSize)
{
  void *ptr;
  if (iSize <= 0) {
     return NULL;
  }
  ptr = (void *)malloc(iSize);
  if(ptr == NULL) {
    printf("Could not allocate memory.\n");
    exit(0);
  }
  return ptr;
}

static inline void* comdCalloc(size_t num, size_t iSize)
{
   return calloc(num, iSize);
}

inline void* comdRealloc(void* ptr, size_t iSize)
{
  ptr = (void *)realloc(ptr, iSize);
  if(ptr == NULL) {
    printf("Could not allocate memory.\n");
    exit(0);
  }
  return ptr;
}

static inline void comdFree(void *ptr)
{
   free(ptr);
}
#endif
