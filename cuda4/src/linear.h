#ifndef __linear_h__
#define __linear_h__

#include <cstdint>

#include "cuda_defines.h"


namespace hw
{
    //mappings for 2d array to simulate fortran arrays in a linear iterations
    //from 0->size you get col and row

    template < int32_t x, int32_t y, int32_t width, int32_t height, int32_t step_x = 1, int32_t step_y = 1 >
    struct linear_2d
    {
        __host__ __device__ static inline int32_t size_x()
        {
            return ((width - x + 1) + step_x - 1) / step_x;
        }

        __host__ __device__ static inline int32_t size_y()
        {
            return ((height - y + 1) + step_y - 1) / step_y;
        }

        __host__ __device__ static inline int32_t size()
        {
            return size_x() * size_y();
        }

        //given index in [0;size] // returns 2d index
        __host__ __device__ static inline int32_t row(uint32_t index)
        {
            return step_y * (index / size_x()) + y;
        }

        __host__ __device__ static inline int32_t col(uint32_t index)
        {
            return step_x * (index % size_x()) + x;
        }
    };

    //mappings for 2d array to simulate fortran arrays in a linear iterations
    //from 0->size you get col and row
    template < int32_t x, int32_t width, int32_t step_x = 1 >
    struct linear_1d
    {
        __host__ __device__ static inline int32_t size_x()
        {
            return ((width - x + 1) + step_x - 1) / step_x;
        }

        __host__ __device__ static inline int32_t size()
        {
            return size_x();
        }

        __host__ __device__ static inline int32_t col(uint32_t index)
        {
            return step_x * (index % size_x()) + x;
        }
    };
}


#endif
