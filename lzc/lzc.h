#ifndef LZC_H_INCLUDED
#define LZC_H_INCLUDED 1

#include <stdio.h>

#if defined(_MSC_VER)
#define LZCEXPORT __declspec(dllexport) 
#else
#define LZCEXPORT
#endif

LZCEXPORT int lz_complexity(int *s, int N);
LZCEXPORT int lz_complexity2(int* s, int N, int threshold);

#endif