/*
This function will calculate the lzc of a sequence probably more quickly than in python
    INPUT:
      *s: int array
      N: int, length of int array
*/
#include "lzc.h"


LZCEXPORT int lz_complexity(int *s, int N)
{
  /* Need to rebuild array?
  for (int i = 0; i < N; i++)
  {
    int seq = s[i]
  }
  */
  int i = 0, k = 1, l = 1;
  int k_max = 1;
  int n = N-1;
  int lzc = 1;

  while(1)
  {
    if (s[i+k-1] == s[l+k-1])
    {
      k += 1;
      if (l+k >= n-1)
      {
        lzc += 1;
        break;
      }
    }
    else
    {
      if (k > k_max)
      {
        k_max = k;
      }
      i += 1;
      if (i == l)
      {
        lzc += 1;
        l += k_max;
        if (l+1 > n)
        {
          break;
        }
        else
        {
          i = 0;
          k = 1;
          k_max = 1;
        }
      }
      else
      {
        k = 1;
      }
    }
  }
  return lzc;
}

// Dummy definition to satisfy Microsoft compiler when using the python extension mechanism to build.
void PyInit_lzc(void)
{

}
