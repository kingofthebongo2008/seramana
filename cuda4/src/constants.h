#ifndef __constants_h__
#define __constants_h__

#include <cstdint>

#include "cuda_defines.h"

namespace hw
{ 
    namespace constants1
    {
        static const int m = 10000;
        static const int n = 100;

        static const float pi     = 3.1415926585f;
        static const float pi2inv = .5f / pi;
    }

    namespace constant_functions
    {
        __host__ __device__ static inline int m()
        {
            return 10000;
        }

        __host__ __device__ static inline int n()
        {
            return 100;
        }

        __host__ __device__ static inline float pi()
        {
            return 3.1415926585f;
        }

        __host__ __device__ static inline float pi2inv()
        {
            return 0.5f / pi();
        }


    }
}




#endif